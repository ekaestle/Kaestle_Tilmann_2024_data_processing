#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 15:32:39 2021

@author: emanuel

most parts of this script are copied from Laura Ermert (ants_2) https://github.com/lermert/ants_2
"""
from mpi4py import MPI
import numpy as np
import os, glob
import sys
import time

from obspy import UTCDateTime, read, Stream, read_inventory, Inventory
from obspy.clients.fdsn.client import Client
from obspy.core.event.catalog import Catalog
from obspy.geodetics import gps2dist_azimuth
from obspy.taup import TauPyModel
from scipy.signal import cheb2ord, cheby2, zpk2sos, sosfilt

import matplotlib.pyplot as plt



comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
print("Hello from rank %g" %rank)
print("MPI size is %g" %size)


"""
PREPROCESSING PARAMETERS
"""
input_directory = "/media/emanuel/Toshiba2TB/SwathD"

input_fileformat = [] # you can specify file endings ["MSEED","sac"], otherwise leave empty []

inventory_directory = "/media/emanuel/Toshiba2TB/SwathD/station_inventory" # xml inventory files

output_directory = "/media/HD3/database_alparray"

overwrite_existing = True

delete_processed_inputfiles = False # delete original files if processed successfully

# remove global events from gcmt catalog using Ekstroems rule
# (includes surface waves circulating several times around the earth)
gcmt_exclude =  False
gcmt_begin   =  '2017,01,01'
gcmt_end     =  '2020,01,01'
minmag       =  6

# removing local events contained in the standard IRIS catalog
# this removes everything from the first arriving phase to the slowest surface wave
event_exclude_cat=True
event_exclude_cat_begin = '2017,01,01' # startime for local catalog search
event_exclude_cat_end =   '2020,01,01' # endtime for local catalog search
event_exclude_cat_minmag = 2.0
event_exclude_cat_lat = 45. # central lat for local catalog search
event_exclude_cat_lon = 12.0 # central lon for local catalog search
event_exclude_cat_radius = 15. #radius in degree

# remove also global events with the same method as for the local catalog
# between the first arriving phase and a signal arriving with 1 km/s, 
# everything is removed
event_exclude_global = True
event_exclude_global_minmag = 7.


min_trace_length = 3*60*60 # minimum length of a trace in seconds
maxgap = 1 # maxiumum allowed gap when merging traces in seconds

# taper, filter and downsample data
downsample = True
sampling_rate = 2.5 # desired sampling rate in Hz, Nyquist freq is half the sampling rate
# freqmax is automatically chosen by the cheby lowpass filter prior to downsampling
freqmin = 1./200 # remove very low frequency signals that are not interesting

# do an instrument correction
instr_correction = True
instr_correction_unit = 'VEL'
instr_correction_prefilt = (1/300, freqmin, 1.25, 2.5) # (f1, f2, f3, f4) or None
# if prefilt is None, it is recommended to set the waterlevel to a lower value (30)
instr_correction_waterlevel = 60.

# remove high energy windows (events and other bursts)
high_e_exclude  =  True
high_e_exclude_level  = 4.0
high_e_exclude_n  = 4
high_e_exclude_std = 4.0
high_e_exclude_winsec = [14400,3600,1800]

"""
END OF PARAMETERS
"""
#%%
"""
PROCESSING FUNCTIONS
"""

def get_event_filter(catalogue):
    """
    Create a time-domain filter removing all events with Mw > 5.6 
    according to GCMT catalogue and 
    the empirical rule of Ekstroem (2001):
    T = 2.5 + 40*(Mw-5.6) [hours]

    catalogue: obspy catalogue object
    """

    event_filter_list = []

    for cat in catalogue[::-1]:
        
        # get origin time
        t_o = cat.origins[0].time
        #print('Original catalog onset: '+t_o.strftime("%Y.%j.%H:%M:%S"))

        t_start = t_o-10

        if len(event_filter_list) > 0:
            if t_o < event_filter_list[-1][1]:
                t_start = event_filter_list[-1][0]
                event_filter_list.pop()
            


        #print('Selected onset: '+t_start.strftime("%Y.%j.%H:%M:%S"))
        # get magnitude
        m = cat.magnitudes[0].mag
        if cat.magnitudes[0].magnitude_type != 'MW':
            raise ValueError('Magnitude must be moment magnitude.')
        if m < 5.6:
            raise ValueError('Event Mw < 5.6: Not ok with Ekstroem event exclusion rule.')

        # determine T
        T = 2.5 * 3600 + 61.8 * (m-5.6) * 3600 

        event_filter_list.append([t_start,t_start+T+10])


    return event_filter_list



def prepare_stream(stream, min_trace_length, maxgap = 1):
    
    for tr in stream:
        if np.isinf(tr.data).any() or np.isnan(tr.data).any() or np.std(tr.data)==0.:
            stream.remove(tr)
    
    # merge overlapping traces that match
    stream._cleanup()    
    
    # remove zero length traces
    for tr in stream:
        if tr.stats.npts <= 0:
            stream.remove(tr)
            
    if len(stream) == 0:
        raise Exception
    
    # check that the sampling rate is the same for all traces
    sampling_rate = tr.stats.sampling_rate
    
    for tr in stream:
        if tr.stats.sampling_rate != sampling_rate:
            raise Exception
        
    # check for overlapping traces that do not match. slice traces
    stream.sort()
    st_temp = Stream()
    st_cleaned = Stream()
    endtime = UTCDateTime('1970,01,01')
    for i,tr in enumerate(stream):
                
        if i > 0:
            endtime = st_temp[0].stats.endtime
            
            # check overlaps
            if endtime > tr.stats.starttime:
                if endtime > tr.stats.endtime:
                    continue
                else:
                    tr = tr.slice(endtime,tr.stats.endtime,nearest_sample=False)
            
            # check gaps
            #if the gap is very small (less than 1 second)
            if tr.stats.starttime - endtime < maxgap:
                if np.abs(tr.data[0] - st_temp[0].data[-1]) < 0.1*np.std(tr.data):
                    st_temp += tr
                    st_temp.merge(method=1,fill_value='interpolate')
                    tr = st_temp[0]
                    st_cleaned.remove(st_cleaned[-1])
        
        st_cleaned += tr
        
        st_temp = Stream(tr)
        
        
    for i in range(len(st_cleaned)-1):
        if st_cleaned[i+1].stats.starttime < st_cleaned[i].stats.endtime:
            print("This should not be possible!")
            print(st_cleaned[0].stats.station,st_cleaned[0].stats.channel,st_cleaned[0].stats.starttime.julday)
            raise Exception

    # remove zero length traces
    for tr in st_cleaned:
        if tr.stats.npts <= 0 or tr.stats.endtime-tr.stats.starttime < min_trace_length:
            st_cleaned.remove(tr)
            
    # slice traces so that they start at a full second if possible
    st_cut = Stream()
    for tr in st_cleaned:
        starttime = np.ceil(abs(tr.stats.starttime))
        st_cut += tr.slice(starttime=UTCDateTime(starttime))
            
    #plt.figure()    
    #for tr in st_cleaned:
    #    tax = np.linspace(abs(tr.stats.starttime),abs(tr.stats.endtime),tr.stats.npts)
    #    plt.plot(tax,tr.data)
    #plt.show()
            
    if len(st_cut) == 0:
        raise Exception

    return st_cut
        


def get_antialias(Fs,freq,maxorder=12):
    # From obspy / L. Ermert's ants_2, maxorder set to 12 (12 in obspy, 8 in 
    # Laura's scripts) and attenuation set to 80 (96 in obspy)
    nyquist = Fs * 0.5
    # rp - maximum ripple of passband, rs - attenuation of stopband
    rp, rs, order = 1, 80, 1e99
    ws = freq / nyquist  # stop band frequency
    wp = ws  # pass band frequency
    # raise for some bad scenarios
    if ws > 1:
        ws = 1.0
        print("** Selected corner frequency is above Nyquist. " + \
              "** Setting Nyquist as high corner.")
    while True:
        if order <= maxorder:
            break
        wp = wp * 0.99
        order, wn = cheb2ord(wp, ws, rp, rs, analog=0)

    z, p, k = cheby2(order, rs, wn,
                     btype='low', analog=0, output='zpk')

    return zpk2sos(z, p, k)
 
    # from scipy import signal
    # b, a = signal.iirfilter(12, [2*np.pi*0.005, 2*np.pi*1.25], rs=80,
    #                         btype='bandpass', analog=True, ftype='cheby2')
    
    # w, h = signal.freqs(b, a, 1000)
    # fig = plt.figure()
    # ax = fig.add_subplot(1, 1, 1)
    # ax.semilogx(w / (2*np.pi), 20 * np.log10(np.maximum(abs(h), 1e-5)))
    # ax.set_title('Chebyshev Type II bandpass frequency response')
    # ax.set_xlabel('Frequency [Hz]')
    # ax.set_ylabel('Amplitude [dB]')
    # ax.axis((0.001, 3, -100, 10))
    # ax.grid(which='both', axis='both')
    # plt.show()



def downsample_stream(stream, freqmin, sampling_rate,
                      taper=True, taper_fraction=0.05):

    stream.detrend('simple') # I think this is better than demean at this point 
    stream.taper(taper_fraction,max_length=int(1.5/freqmin))
    
    freqmax = 0.5*sampling_rate
    
    stream.filter("bandpass", freqmin=freqmin, freqmax=freqmax, corners=3, zerophase=True)
    # butterworth filter tapers off too smoothly to make sure there is no
    # aliasing effect

    anti_alias = get_antialias(stream[0].stats.sampling_rate,freqmax)
    for tr in stream:
        # zerophase cheby lowpass from L. Ermert's ants_2
        firstpass = sosfilt(anti_alias, tr.data)
        tr.data = sosfilt(anti_alias, firstpass[::-1])[::-1]

    for tr in stream:
        
        # sometimes the sampling rate is something like 20.0000115
        # in this case, decimation will fail unless the sampling rate is set
        # to 20.0
        if tr.stats.sampling_rate > 1:
            # if the trace with the rounded sampling rate is shifted by less
            # than one sample from the original sampling rate
            if np.abs(tr.stats.npts/tr.stats.sampling_rate - tr.stats.npts/np.round(tr.stats.sampling_rate,2)) < 1./sampling_rate:
                tr.stats.sampling_rate = np.round(tr.stats.sampling_rate,2)
        
        if tr.stats.sampling_rate <= sampling_rate:
            continue
        
        downsample_interval = tr.stats.sampling_rate/sampling_rate
        
        # make sure the downsampling interval is an integer
        if np.abs(np.round(downsample_interval) - downsample_interval) > 0.:
            print("downsampling rate is badly chosen. trace sampling rate:",tr.stats.sampling_rate)
            raise Exception
        
        tr.decimate(int(downsample_interval),no_filter=False)
        
        
    stream.detrend(type='linear')
    
    return stream



def remove_events(stream, event_filter, event_catalog, inv, 
                  min_trace_length, ofid, minmag=5.6):
    

    if event_filter is not None:    

        for tr in stream:
            tr.detrend('demean')
        
        t_total = 0.0
        for trace in stream:
            t_total += trace.stats.npts
        
        for quake_window in event_filter:
            stream.cutout(starttime=quake_window[0],
                          endtime=quake_window[1])
        
        t_kept = 0.0
        for trace in stream:
            t_kept += trace.stats.npts
        
        print('* Excluded all events in GCMT catalogue with Mw >=' +
              str(minmag),
              file=ofid)
        print('* Lost %g percent of original traces'
              % ((t_total - t_kept) / t_total * 100), file=ofid)

        
    if event_catalog is not None:
        
        model = TauPyModel(model="iasp91")
        
        for tr in stream:
            tr.detrend('demean')

        t_total = 0.0
        for trace in stream:
            t_total += trace.stats.npts
             
        for event in event_catalog:
            # get origin time
            t0 = event.origins[0].time
            lon0 = event.origins[0].longitude
            lat0 = event.origins[0].latitude
            depth0 = event.origins[0].depth/1000.
            coords = inv.get_coordinates(stream[0].id,datetime=st[0].stats.starttime)
            data_start = stream[0].stats.starttime
            if t0 < data_start-24*60*60.:
                continue
            data_end = stream[-1].stats.endtime
            if t0 > data_end:
                continue
            dist = gps2dist_azimuth(lat0,lon0,coords["latitude"],
                                    coords["longitude"])[0]/1000.
            # p_arrival is ordered from fastest to slowest arrival
            p_arrival = model.get_travel_times(source_depth_in_km=depth0,
                                               distance_in_degree=dist/111.19,
                                               phase_list=["P","Pdiff","PKP"])
            if len(p_arrival)==0:
                tcut1 = t0
            else:
                tcut1 = t0 + p_arrival[0].time - 60.0 #60s before first arrival
            if tcut1<t0:
                tcut1 = t0
            tcut2 = t0 + dist/1.0 + np.max([dist/1.0*0.1,60]) #surface-wave arrival of 1km/s plus 10% (min 60s)
            stream.cutout(starttime=tcut1,endtime=tcut2)


        t_kept = 0.0
        for trace in stream:
            t_kept += trace.stats.npts
        
        print('* Excluded all events in event catalogue.', file=ofid)
        print('* Lost %g percent of original traces' %((t_total-t_kept)/t_total*100), file=ofid)


    # remove zero length traces
    for tr in stream:
        if tr.stats.npts <= 0 or tr.stats.endtime-tr.stats.starttime < min_trace_length:
            stream.remove(tr)
            
    if len(stream) == 0:
        print('** No traces left after event removal.',file=ofid)
        raise Exception
    



def remove_high_energy_windows(stream,windows,n_compare,min_trace_length,\
                               factor_enrg=1.,taper_perc=0.05,thresh_stdv=1.,
                               ofid=None):

    """
    Modified from L. Ermerts ants_2
    windows = high_e_exclude_winsec
    n_compare = high_e_exclude_n
    factor_enrg = 5
    thresh_stdv = 5
    taper_perc=0.05
    A really somewhat complicated way to try and get rid of earthquakes and other high-energy bursts. 
    A sort of coarse multiwindow-trigger; I haven't found a better (comparatively fast) way so far.
    High-energy events will be cut out.
    Operates directly on the trace.
    """
    
    print('* Removing high energy windows ',file=ofid)
    cutout_windows = []
    pctg_cut = 0.
    
    # should I run this twice as in the original version from Laura?
    # I think it makes no sense when I use the cutout function as I do here
    for k in range(1):
        
        for trace in stream:
        
            if trace.stats.npts == 0:
                print("zero length trace is not allowed at this stage")
    
            windows.sort() # make sure ascending order
            cutout_times = 0.
            
            testtrace = trace.copy()
        
            for proc_steps in trace.stats.processing:
                if "bandpass" in proc_steps:
                    freq_min = float(proc_steps.split("'freqmin': ")[-1].split(",")[0])
                    break
            else:
                freq_min = 0.01
                freq_max = 1.0
                testtrace.taper(type='cosine',max_percentage = taper_perc,
                                max_length=int(1.5/freq_min))
                testtrace.filter('bandpass',freqmin = freq_min,freqmax = freq_max,
                                  corners = 3, zerophase = True)
            
            for iwin,win in enumerate(windows): 
            # initialize arrays of subtrace values (energy, standard deviation) for each window length; maybe just use an expandable list (max. a couple of hundred values)
                enrg = []
                stdv = []
                t = []
                window_times = []
            # fill those arrays
                t0 = trace.stats.starttime
                last_window = False
                while True:
                    
                    subtr = testtrace.slice(starttime=t0,endtime=t0+win-1).data
                    
                    enrg.append(np.sum(np.power(subtr,2))/win)
                    subwin = int(win/3)
                    [a,b,c] = [ np.std(subtr[0:subwin]),
                                np.std(subtr[subwin:2*subwin]),
                                np.std(subtr[2*subwin:3*subwin]) ]  
                                      
                       
                    stdv.append(np.max([a,b,c])/(np.min([a,b,c])+sys.float_info.epsilon))
                    
                    window_times.append([t0,t0+win-1])
                                        
                    t.append((t0+win/2).strftime('%s'))
                    t0 += int(win*0.7) #30% overlap
                    
                    if last_window:
                        break                
                    if t0 > trace.stats.endtime-win:
                        t0 = trace.stats.endtime-win
                        last_window = True

                # count how many windows are excluded on counter
                # step through the array enrg, stdv; this should be relatively fast as array are not particularly long
        
                         
                for i in range(len(enrg)):
  
                    mean_enrg = np.mean(np.roll(enrg,n_compare-i)[:n_compare])
                    
                    # i0 = i - n_compare if i>=n_compare else 0
                    # i1 = i0 + n_compare if i0+n_compare<=len(enrg) else len(enrg)
                    # if i0+1>=i1:
                    #     continue
                    # mean_enrg = np.mean(enrg[i0:i1])
                                        
                    if ((enrg[i] > factor_enrg * mean_enrg and stdv[i] > thresh_stdv) or
                        (enrg[i] > 10 * mean_enrg)):
        
                        t0 = window_times[i][0]
                        t1 = window_times[i][1]

                        cutout_windows.append([t0,t1])
        
                        cutout_times += t1-t0
   
                        
            # percentage of data that was cut:
            pctg_cut += cutout_times / (trace.stats.endtime-trace.stats.starttime) * 100
            
            #if k == 0:
            #    print("first round, %.1f percent cut" %(float(np.sum(weight==0.))/float(trace.stats.npts)*100))
            #else:
            #    print("second round, %.1f percent cut" %(float(np.sum(weight==0.))/float(trace.stats.npts)*100))
                    
            
        for cut_win in cutout_windows:
            stream.cutout(cut_win[0],cut_win[1])
                

    # remove zero length traces
    for tr in stream:
        if (tr.stats.npts <= 0 or
            tr.stats.endtime - tr.stats.starttime < min_trace_length):
            stream.remove(tr)
            
    # display a summary of how much was kept 
    if True:# pctg_cut > 0:
        print('* cut %d percent of data from trace: ' %(pctg_cut),file=ofid)
        
    if len(stream) == 0:
        print('** No traces left after removing high energy windows. ',file=ofid)
        raise Exception
        



def remove_response(stream, inventory, pre_filt, waterlevel, unit, ofid, 
                    taper_perc = 0.05):   

    if isinstance(inventory, dict):
        print("currently only stationxml is supported.")
        raise Exception()

    elif isinstance(inventory, Inventory):

        freq_min = 0.00001
        if pre_filt is not None:
            freq_min = pre_filt[1]
        else:
            for proc_steps in stream[0].stats.processing:
                if "bandpass" in proc_steps:
                    freq_min = float(proc_steps.split("'freqmin': ")[-1].split(",")[0])
                    break
            else:
                print("Warning! It seems no pre-filtering and no bandpass was "+
                      "applied prior to removing the instrument response. " +
                      "may cause unwanted effects.")

        stream.detrend('demean')

        stream.taper(type='cosine', max_percentage = taper_perc,
                     max_length=int(1.5/freq_min))

        stream.remove_response(inventory=inventory,
                               pre_filt=pre_filt,
                               water_level=waterlevel,
                               taper = False,
                               output=unit)
    else:
        msg = 'No inventory or seedresp found.'
        raise ValueError(msg)
        
        
            
"""
END PROCESSING FUNCTIONS
"""            
            


#%%
"""

This script preprocesses the MSEED files in the input directories 
specified in the input file.
 

"""
# Create output directory, if necessary

 
if rank == 0 and not os.path.exists(output_directory):
    os.makedirs(output_directory)
if rank == 0:
    print("Starting preprocessing")

comm.Barrier()

event_filter = None
if gcmt_exclude:

    if rank == 0:
        print("getting GCMT catalog...")
        c = Client()
        cata = c.get_events(starttime=UTCDateTime(gcmt_begin),
            endtime=UTCDateTime(gcmt_end),catalog='GCMT',
            minmagnitude=minmag)

        event_filter = get_event_filter(cata)

    # communicate event_filter (would it be better 
    # if every rank sets it up individually?)
    event_filter = comm.bcast(event_filter,root=0)

event_cat = None
if event_exclude_cat or event_exclude_global:
    
    event_cat = Catalog()
    
    if rank == 0 and event_exclude_cat:
        print("getting event catalog...")
        c = Client()
        event_cat.extend(c.get_events(
                starttime=UTCDateTime(event_exclude_cat_begin),
                endtime=UTCDateTime(event_exclude_cat_end),
                #catalog=catalog,
                minmagnitude=event_exclude_cat_minmag,
                latitude=event_exclude_cat_lat,
                longitude=event_exclude_cat_lon,
                maxradius=event_exclude_cat_radius))
        
    if rank == 0 and event_exclude_global:
        event_cat.extend(c.get_events(starttime=UTCDateTime(event_exclude_cat_begin),
            endtime=UTCDateTime(event_exclude_cat_end),catalog='GCMT',
            minmagnitude=event_exclude_global_minmag))
    
    event_cat = comm.bcast(event_cat,root=0)
    print(len(event_cat),"events in earthquake catalog.")

#%% PROCESSING LOOP

# processing report file
sys.stdout.flush()
output_file = os.path.join(output_directory,
    'processing_report_rank%g.txt' %rank)


if os.path.exists(output_file):
    ofid = open(output_file,'a')
    print('UPDATING, Date:',file=ofid)
    print(time.strftime('%Y.%m.%dT%H:%M'),file=ofid)
else:
    ofid = open(output_file,'w')
    print('PROCESSING, Date:',file=ofid)
    print(time.strftime('%Y.%m.%dT%H:%M'),file=ofid)


while True:
    
    #- Find input files
    file_endings = []
    for fformat in input_fileformat:
        file_endings.append(fformat.split(".")[-1].lower())
    input_fileformat = file_endings
    content = []
    if rank==0:
        content = glob.glob(os.path.join(input_directory,"**/*"),recursive=True)
        # remove folders from the list (everything that is not a file)
        accepted_filepaths = []
        for filepath in content:
            if os.path.isfile(filepath) and not (filepath.endswith("xml")):
                accepted_filepaths.append(filepath)
        content = accepted_filepaths
        # remove files that do not fit the input_fileformat pattern
        if len(input_fileformat) > 0:
            accepted_filepaths = []
            for filepath in content:
                if content.split(".")[-1].lower() in input_fileformat:
                    accepted_filepaths.append(filepath)
            print(len(content)-len(accepted_filepaths),"files are being ignored because of the input_fileformat parameter.")
            content = accepted_filepaths
        print(len(content), "files found to be processed")
        np.random.shuffle(content)
    content = comm.bcast(content,root=0)

    
    # find inventory files
    inventory_filelist = []
    if rank == 0:
        inventory_filelist = glob.glob(os.path.join(inventory_directory,"**/*.xml"),recursive=True)
    inventory_filelist = comm.bcast(inventory_filelist,root=0)
    
    # select input files for this rank
    content = np.sort(content)
    content = content[rank::size]
    
    # Loop over input files
    t0 = time.time()
    for i,filepath in enumerate(content):

        #if not ("ABTA" in filepath or "BIOA" in filepath):
        #    continue
        #print(filepath)

        #if not ("DAVOX" in filepath and "BHN.D.2018.114" in filepath):# and not ("BRMO" in filepath):
        #    continue
        #pause
        #if not ("HZ" in filepath):
        #    continue
        #if "D013" in filepath and "HHN" in filepath and "334" in filepath:
        #    pause
        #else:
        #    continue
            
        if i%1000==0:
            print("MPI Rank",rank,"    ",i,"/",len(content))
            if i>0:
                print("time left: %.1f mins" %((len(content)-i)/(i+1)*(time.time()-t0)/60.))
        
        print('-------------------------------------',file=ofid)
        print('Attempting to process:',file=ofid)
        print(os.path.basename(filepath),file=ofid)
        
        outfile = filepath.replace(input_directory,output_directory)
          
        if os.path.isfile(outfile):
            if not overwrite_existing:
                print('** File already exists (skipping):',outfile)
                print('** File already exists, skipping:',file=ofid)
                print('** %s' %filepath,file=ofid)
                continue
            else:
                print('** File already exists (overwriting):',outfile)
                print('** File already exists, overwriting.',file=ofid)
        
        outdir = os.path.dirname(outfile)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        
        try:
            st = read(filepath)
        except:
             print('** Problem opening file, skipping: ',file=ofid)
             print('** %s' %filepath,file=ofid)
             continue
    
        if len(st) == 0:
            print('** No data in file, skipping: ',file=ofid)
            print('** %s' %filepath,file=ofid)
            continue
        
        st_original = st.copy()    
        # get the inventory file
        station = st[0].stats.station
        network = st[0].stats.network
        try:
            for tr in st:
                if tr.stats.network != network or tr.stats.station != station:
                    print('** File contains more than one station data, skipping: ',file=ofid)
                    print('** %s' %filepath,file=ofid)
                    raise Exception()
            inv_filepaths = []
            for inv_filepath in inventory_filelist:
                if (inv_filepath.split("/")[-1].split(".")[0] == network and
                    inv_filepath.split("/")[-1].split(".")[1] == station):
                    inv_filepaths.append(inv_filepath)
            if len(inv_filepaths) == 1:
                inv_filepath = inv_filepaths[0]
            else:
                raise Exception("more than one inventory file for station",network,station,inv_filepaths)
            inventory = read_inventory(inv_filepath)
        except:
            print('** Could not find station inventory, skipping: ',file=ofid)
            print('** %s' %filepath,file=ofid)
            continue
            
        try:
            st = prepare_stream(st, min_trace_length, maxgap=maxgap)
        except:
            print('** Problems preparing stream: ',file=ofid)
            print('** %s' %filepath,file=ofid)
            continue
    
    
        try:
            remove_events(st, event_filter, event_cat, inventory,
                          min_trace_length, ofid, minmag=minmag)
        except:
            print('** Problems removing earthquake events from stream: ',file=ofid)
            print('** %s' %filepath,file=ofid)
            continue
    

        if downsample:         
            #try:
            st = downsample_stream(st, freqmin, sampling_rate,
                                   taper=True, taper_fraction=0.05)
            #except:
            #    print('** Problems downsampling stream: ',file=ofid)
            #    print('** %s' %filepath,file=ofid)
            #    continue
        
        st_after_eventexclude = st.copy()
        if high_e_exclude:
            try:
                remove_high_energy_windows(st,high_e_exclude_winsec,high_e_exclude_n,min_trace_length,
                                           factor_enrg=high_e_exclude_level,
                                           taper_perc=0.05,
                                           thresh_stdv=high_e_exclude_std,ofid=ofid)
            except:
                print('** Problems excluding high energy windows stream: ',file=ofid)
                print('** %s' %filepath,file=ofid)
                continue
        
        
        # traces are now filtered (from downsampling)
        # traces were tapered while downsampling but not tapered after event
        # exclude. remove response will taper the data before the response is
        # removed.
        
        if instr_correction:
            try:
                remove_response(st, inventory, instr_correction_prefilt,
                                instr_correction_waterlevel,
                                instr_correction_unit, ofid)
            except:
                print('** Problems removing response: ',file=ofid)
                print('** %s' %filepath,file=ofid)
                continue
        
        invalid = False
        for tr in st:
            if np.isinf(tr.data).any() or np.isnan(tr.data).any() or np.std(tr.data)==0.:
                print('** invalid values after processing: ',file=ofid)
                print('** %s' %filepath, file=ofid)
                invalid = True
        if invalid:
            continue
                
        
        
        try:
            for tr in st:
                tr.data = tr.data.astype('float32')
                tr.stats.mseed.encoding = 'FLOAT32'
            st.write(outfile,format='MSEED')
            print('** success ',file=ofid)
            #st.write(outfile,format="PICKLE")
        except:
            print('** Problems writing stream: ',file=ofid)
            print('** %s' %filepath,file=ofid)
            continue
    
    
        if i%1000==0:
            tstart = UTCDateTime(year=(st[0].stats.starttime+60).year,
                                 julday=(st[0].stats.starttime+60).julday)
            plt.ioff()
            fig = plt.figure(figsize=(16,12))
            plt.subplot(311)
            for tr in st_original:
                timeax = np.linspace(tr.stats.starttime-tstart,
                                     tr.stats.endtime-tstart,tr.stats.npts)
                plt.plot(timeax,tr.data,'b')
                plt.xlim(0,24*60*60)
            plt.subplot(312)
            for tr in st_after_eventexclude:
                timeax = np.linspace(tr.stats.starttime-tstart,
                                     tr.stats.endtime-tstart,tr.stats.npts)
                plt.plot(timeax,tr.data,'b')
                plt.xlim(0,24*60*60)
            plt.subplot(313)
            for tr in st:
                timeax = np.linspace(tr.stats.starttime-tstart,
                                     tr.stats.endtime-tstart,tr.stats.npts)
                plt.plot(timeax,tr.data,'b')     
                plt.xlim(0,24*60*60)
            plt.savefig(os.path.join(output_directory,outfile.split("/")[-1]+"_2.jpg"),dpi=200,bbox_inches='tight')
            plt.close(fig)   
        
        if delete_processed_inputfiles:
            os.remove(filepath)
    
        ofid.flush()
        
        time.sleep(0.001)
    
    break
    print("Rank %g has completed processing, waiting eight hours and checking for new files..." %rank)
    print("Sleeping until",time.ctime(time.time()+8*60*60))
    time.sleep(8*60*60)
    
 
ofid.close()

print("Rank %g has completed processing." 
    %rank,file=None)
