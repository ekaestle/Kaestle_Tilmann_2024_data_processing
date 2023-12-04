# -*- coding: utf-8 -*-
"""
Created in July 2015

@author: emanuel

This script will process all previously created cross correlation files in 'spectra_path'.
It is assumed that the naming scheme for the CC files is according to the create_cc.py script.
Otherwise please adapt command below, so that station names are recognized correctly from
the file name and can be found in the station list file.
If you encounter any difficulties, do not hesitate to contact me: Emanuel Kaestle @ FU Berlin
"""


from mpi4py import MPI
mpi_comm = MPI.COMM_WORLD
mpi_rank = mpi_comm.Get_rank()
mpi_size = mpi_comm.Get_size()

import glob,os,sys
import numpy as np
import noise
import datetime
import random
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy import optimize
#from obspy.geodetics import gps2dist_azimuth
from obspy.signal.filter import bandpass
#from obspy import read
import pickle

""" you can run this program also by calling
python process_spectra.py stat1 stat2
it will then only process stat1 and stat2 and create several plots """

""" USER PARAMETERS"""
max_period = 999. # don't use data above that period
min_period = 1 # don't use data below that period
min_wavelength_criterion = 1 # don't pick data at less than one wavelength inter-station distance
max_wavelength_criterion = 50 # don't pick data at more than x wavelengths inter-station distance
min_vel=1.5 # minimum velocity for the desired wave type (Rayleigh/Love)
max_vel=5.5 # maximum velocity ...
min_common_days = 180 # recommended to be at least half a year
#min_amplitude = 0.05 # minimum allowed amplitude of the cc spectrum with respect to its maximum (should be tested, otherwise set to 0)
pick_thresh = 1.8
horizontal_polarization=False # If True, TT or RR cross-correlations are calculated (additional J2 term). If False, only ZZ
spectra_path = "/media/HD3/cross_correlation_spectra_3hrs/ZZ/" # path where spectra files are stored (filenames are expected to have both station names, as the ones created by create_ccs.py. Otherwise adapt stat1 = fname.split("/")[-1].split("_")[0] below.)
phase_vel_path='phasevel_curves_monthly_ZZ' # path where to save the cross_correlation spectra
ref_curve = np.loadtxt("../rayleigh/Average_phase_velocity_rayleigh") # reference curve for picking
velocity_filter = True # do a velocity filter between min and max velocity, recommended
#check_pos_neg_lagtime_correlations = True # accept dispersion curves only if the ones picked from the positive and negative lag time cross correlations are similar
plotresult = False
# SECTION ABOUNT SNR FILTERS. This is a filter applied only at the very end, discarding parts of the phase-velocity curves where the SNR is lower than a certain threshold
# this has not been tested a lot yet. I recommend to keep the filters on False for most applications.
# snr filters are only active in combination with the velocity filter
# snr filter in the time domain:
snr_filter_time_domain = False
snr_filter_threshold_td = 2
#snr filter in the freq domain
snr_filter_frequency_domain = False
snr_filter_threshold_fd = 1 # see Sadeghisorkhani 2017 (computers and geosciences), appendix; only important when snr_filter_frequency_domain=True
# used for both freq and time domain snr filters. Checks the SNR in certain frequency intervals.
# define intervals where to check the snr ratio, e.g.: 1./np.array([200,100,50,35,20,10,6,3])
intervals = np.append(0,np.logspace(np.log10(0.01),np.log10(0.3),10))
"""END OF USER DEFINED PARAMETERS """
counter = 0
countermax = 999999999999999
""" ADDITIONAL PARAMETERS AND NOTES
- It may be interesting to modify also other parameters that influence the result
of the picking procedure. For example the filter width and height for the 
get_smooth_pv function.
- If the velocity filter is switched off, smoothing will be applied to the
spectrum before picking. If smoothing should be switched off entirely, please
modify values in the picking funtion below.
- You can try also the classical zero-crossing picking procedure. In order to 
do so, please replace get_smooth_pv with extract_phase_velocity. In this case,
remove also the internal try/except loop and remove input values filt_width, 
filt_height, x_step, pick_threshold. See documentation of noise.py.
- the file name recognition pattern may have to be changed if your cross
correlation files have a different naming scheme. Please check if you encounter
any errors.
- feel free to contact me if you have questions: emanuel.kaestle@fu-berlin.de
"""

if (velocity_filter==False and snr_filter_time_domain) or (velocity_filter==False and snr_filter_time_domain):
    print("snr filter only works with velocity filter activated!")
    sys.exit()

def running_mean(x, N):
    if N%2 == 0:
        N+=1
    x = np.insert(x,0,np.ones(np.int(N/2))*x[0])
    x = np.append(x,np.ones(np.int(N/2))*x[-1])
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / N 
    
def snr_cc(symmetric_cc,df,distance,cmin,cmax,intervals,plotting=False):
    snr_filt = np.zeros(len(intervals)-1)
    if plotting:
        plt.figure()
        plt.subplot(len(intervals)+1,1,1)
        plt.plot(symmetric_cc)
        ax = plt.gca()
        ax.set_yticklabels([])
    for i in range(len(intervals)-1):
        lim1 = intervals[i]
        lim2 = intervals[i+1]
        signal = bandpass(symmetric_cc,lim1,lim2,df)
        idx1 = int(dist/cmax*df)
        idx2 = int(dist/cmin*df)
        snr_filt[i] = np.mean(np.abs(signal[idx1:idx2])) / np.mean(np.abs(np.append(signal[0:idx1],signal[idx2:int(len(signal)/2)])))
        if plotting:
            plt.subplot(len(intervals)+1,1,i+2)
            plt.plot(signal)
            plt.plot([idx1,idx1],[np.min(signal),np.max(signal)],'r')
            plt.plot([idx2,idx2],[np.min(signal),np.max(signal)],'r')
            if lim1==0.:
                lim1 = 1./999
            plt.text(idx2,np.max(signal)/2.,"%d - %ds SNR: %.1f" %(1./lim2,1./lim1,snr_filt[i]))
            ax = plt.gca()
            ax.set_yticklabels([])
    if plotting:
        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
    return snr_filt

        
if np.mean(ref_curve[:,1])>1000:
    ref_curve[:,1]/=1000.
    
statpair = False
try:
    station1=sys.argv[1]
    station2=sys.argv[2]
    statpair = True
    plotresult=True
except:
    pass

""" for testing """
station1 = 'CH.DAVOX'
station2 = 'ZS.D043'
statpair = True
plotresult=True
""" """

    
phase_vel_path_lq = os.path.abspath(phase_vel_path)+"_quality_check_failed/"
figure_path = os.path.abspath(phase_vel_path)+"_figures/"
figure_path_failed = os.path.abspath(phase_vel_path)+"_figures_quality_check_failed/"

spectralist = []
if mpi_rank==0:
    if not os.path.exists(phase_vel_path):
        os.makedirs(phase_vel_path)
    if not os.path.exists(phase_vel_path_lq):
        os.makedirs(phase_vel_path_lq)
    if not os.path.exists(figure_path):
        os.makedirs(figure_path)
    if not os.path.exists(figure_path_failed):
        os.makedirs(figure_path_failed)
    
    spectralist = glob.glob(os.path.join(spectra_path,"*")) 
    phase_vel_curves = glob.glob(os.path.join(phase_vel_path,"*"))
    phase_vel_curves += glob.glob(os.path.join(phase_vel_path_lq,"*"))
    print(len(phase_vel_curves),"already calculated phase-velocity curves have been found.")
    
    if not statpair:
        for fpath in phase_vel_curves:
            fname = fpath.split("/")[-1].replace(".curve",".pkl")
            fname_remove = os.path.join(spectra_path,fname)
            try:
                spectralist.remove(fname_remove)
            except:
                print("file appears twice:",fname_remove)
                print("check for the file in pathes:")
                print(phase_vel_path)
                print(phase_vel_path_lq)
                pass
    
    random.shuffle(spectralist)
    
spectralist = mpi_comm.bcast(spectralist,root=0)   
spectralist = spectralist[mpi_rank::mpi_size]

avg_wavelength = interp1d(ref_curve[:,0],ref_curve[:,1]/ref_curve[:,0],bounds_error=False,fill_value='extrapolate')
ref_curve_fu = interp1d(ref_curve[:,0],ref_curve[:,1])
quality_check_collection = []


time_start=datetime.datetime.now()
print("MPI Rank",mpi_rank,"processing",len(spectralist),"cross-correlation spectra.")
m=1
fails=0
for i,fname in enumerate(spectralist):
    if counter>countermax:
        print("STOPPING: Maximum counter number reached")
        break
    #if np.sum(["ZS" in part for part in fname.split("_")])!=2:
    #    continue
    if statpair:
        if not (station1 in fname and station2 in fname):
            continue
        else:
            print(fname)
    if plotresult:
        print("\n # # # # # # # # # # # # # # # # \n",fname)
    if float(i)/len(spectralist) > m*0.01:
        time_passed = (datetime.datetime.now()-time_start).total_seconds()
        time_left = (100.-m)/m * time_passed
        print("MPI Rank",mpi_rank,m*1,"percent processed, estimated time left:",np.round(time_left/60.),"minutes")
        print(fails,"/",i,"cross-correlations could not be processed.")
        m+=1
        
    try:
        with open(fname, "rb") as f:
            ccdict = pickle.load(f)
    except:
        print("\nError: could not read file",fname,"\n\n")
        continue
    if len(ccdict['corrdays']) < min_common_days:
        #print("    ",os.path.basename(fname),"skipping: only %d correlation days" %(len(ccdict['corrdays'])))
        fails += 1
        continue
    spectrum = ccdict['spectrum']
    freqax = ccdict['freq']
    stat1 = ccdict['station1']['id']
    stat2 = ccdict['station2']['id']
    dist = ccdict['dist']
    
    #if dist>100 and not(statpair):
    #    continue
    
    az = ccdict['az']
    baz = ccdict['baz']   
    dt = 1./freqax[-1]/2.
    overlap = float(fname.split("ovlap_")[-1].split(".")[0])
    
    tcorr = np.fft.irfft(spectrum)
    
    if velocity_filter:
        freqax,corr_filtered,timeax,tcorr_filtered = noise.velocity_filter(
            freqax,spectrum,dist,velband=(max_vel+0.1,max_vel,min_vel,0.5*min_vel),
            return_all=True)
        TCORR_full = np.fft.fftshift(tcorr_filtered)
    else:
        timeax = (np.arange(len(tcorr))-len(tcorr)/2.)*dt
        TCORR_full = tcorr
    winlength = np.max(timeax)*2

#    # check whether ray runs through the Ivrea body:
#    s = np.arange(0,dist,20) # 20 km sampling
#    lats = np.linspace(ccdict['station1']['latitude'],ccdict['station2']['latitude'],len(s))
#    lons = np.linspace(ccdict['station1']['longitude'],ccdict['station2']['longitude'],len(s))
#    areacheck = False
#    if ((lats>44.)*(lats<46.)*(lons>6.5)*(lons<8.)).any() and dist<3000.:
#        areacheck = True
#    else:
#        continue
#    if pobasin and not("po_basin" in phase_vel_path):
#        phase_vel_path+="_po_basin"
#        if not os.path.exists(phase_vel_path):
#            os.makedirs(phase_vel_path)
#    elif not(pobasin) and "po_basin" in phase_vel_path:
#        phase_vel_path = phase_vel_path[:-9]


#    dt = tr.stats.sac['delta']
#    freqax = np.fft.rfftfreq(len(crosscorr),dt)
#    id1,id2 = fname.split("/")[-1].split("--")
#    net1,sta1,loc1,cha1 = id1.split(".")
#    stat1 = net1+"."+sta1
#    net2,sta2,loc2,cha2,_ = id2.split(".")
#    stat2 = net2+"."+sta2
#    #stat1 = fname.split("/")[-1].split(".")[0:2].join()
#    #stat2 = fname.split("/")[-1].split(".")[2]
#    #dist = float(fname.split("/")[-1].split("_")[4])/1000.
#    try:
#        dist,az,baz = gps2dist_azimuth(stat_dict[stat1][0],stat_dict[stat1][1],stat_dict[stat2][0],stat_dict[stat2][1])
#    except:
#        print("station information not found in statlist! aborting...")
#        continue
#    dist/=1000.
    #dt = 1./spectrum[-1,0]/2.

    # determine min/max frequency dependent on the wavelength criteria and the
    # average wavelength (avg wavelength from reference curve)
    interp_fn2 = lambda x: avg_wavelength(x)-dist/min_wavelength_criterion
    x0 = ref_curve[np.abs(interp_fn2(ref_curve[:,0])).argmin(),0]
    try:
        f0 = optimize.newton(interp_fn2,x0,tol=1e-4)
    except:
        f0 = x0
    interp_fn2 = lambda x: avg_wavelength(x)-dist/max_wavelength_criterion   
    x0 = ref_curve[np.abs(interp_fn2(ref_curve[:,0])).argmin(),0]
    try:
        f1 = optimize.newton(interp_fn2,x0,tol=1e-4)
    except:
        f1 = x0


    min_freq = np.max([f0,1./max_period])
    max_freq = np.min([f1,1./min_period]) 
        
    #checking data quality: analysing forward, backward and symmetric cross correlations    
    
    FW = TCORR_full[:int(len(TCORR_full)/2)+1]
    FW = np.append(FW,FW[1:-1][::-1])
    
    BW = TCORR_full[int(len(TCORR_full)/2):]
    BW = np.append(BW,TCORR_full[0])
    BW = np.append(BW[::-1],BW[1:-1])
    
    # evaluate (time-domain) forward and backward spectrum separately to check the data quality
    types = ["symmetric","forward","backward"]
    for entry in ccdict:
        if "spectrum." in entry:
            types.append(entry)
           
    speclist = []
    freqlist = []
    for ispec,subspec in enumerate(types):           

        freqax = ccdict['freq']
        if subspec == "forward":
            TCORR = FW
            CORR = np.fft.rfft(TCORR)
        elif subspec == "backward":
            TCORR = BW
            CORR = np.fft.rfft(TCORR)
        elif subspec == "symmetric":
            TCORR = TCORR_full
            CORR = np.fft.rfft(TCORR)
        else:
            nwins = ccdict[subspec.replace("spectrum","no_wins")]
            days_per_month = winlength*nwins*(1-overlap)/60/60/24
            if days_per_month < 20:
                continue
            _s,year,month = subspec.split(".")
            spec_month = ccdict[subspec]
            freqax = np.linspace(0,freqax[-1],len(spec_month))
            
            if velocity_filter:
                _f,CORR,_t,_tcorr = noise.velocity_filter(
                    freqax,spec_month,dist,
                    velband=(max_vel+0.1,max_vel,min_vel,0.5*min_vel),
                    return_all=True)
            else:
                CORR = spec_month
        
        freqlist.append(freqax)
        speclist.append(CORR)
        
        if subspec == 'symmetric':
        
            # no velocity filter applied
            tcorr = np.fft.irfft(spectrum)
            
            idx_tmin = int((dist/max_vel)/dt) 
            idx_tmax = int((dist/min_vel)/dt)
            overall_snr = np.mean(np.mean(np.abs(tcorr[idx_tmin:idx_tmax]))/np.mean(np.abs(np.append(tcorr[0:idx_tmin],tcorr[idx_tmax:int(len(tcorr)/2)]))))
    
            """ applying SNR filter"""
            if snr_filter_time_domain:
                snr_filt_td = snr_cc(tcorr,1./dt,dist,min_vel,max_vel,intervals,plotting=statpair) > snr_filter_threshold_td
            if snr_filter_frequency_domain:
                a_thresh = np.sqrt((1+snr_filter_threshold_fd**2)/((len(spectrum)*2-2)/(0.5*dist*(1./min_vel-1./max_vel))+snr_filter_threshold_fd**2))
                Psignal = np.abs(spectrum)**2
                Psignal = running_mean(Psignal,21)
                Pnoise = np.abs(CORR)**2
                Pnoise = running_mean(Pnoise,21)
                Pnoise[Pnoise>Psignal] = Psignal[Pnoise>Psignal]
                a = np.sqrt(Pnoise/Psignal)
                snr_filt_fd=[]
                for i in range(len(intervals)-1):
                    lim1 = intervals[i]
                    lim2 = intervals[i+1]
                    if lim1==0.:
                        lim1=1./999.
                    if np.sum(a[(lim1<=freqax)*(freqax<lim2)]<a_thresh)>1:
                        snr_filt_fd.append(False)
                    else:
                        snr_filt_fd.append(True)      
                if statpair:
                    plt.figure()
                    plt.plot(freqax,Psignal/np.max(Psignal),label='PSD real cc-spectrum with velocity filter')
                    plt.plot(freqax,Pnoise/np.max(Pnoise),label='PSD real cc-spectrum')
                    plt.plot(freqax,a,label='PSDnoise/PSDsignal')
                    #plt.plot(freqax,running_mean(a,31))
                    plt.plot(freqax,np.ones(len(spectrum))*a_thresh,'--',label='SNR=%.1f threshold' %snr_filter_threshold_fd)
                    for i in range(len(intervals)-1):
                        lim1 = intervals[i]
                        lim2 = intervals[i+1]
                        plt.plot([lim1,lim1],[0,1],'k',lw=0.2)
                        plt.plot([lim2,lim2],[0,1],'k',lw=0.2)
                        if snr_filt_fd[i]:
                            plt.text(lim1,0.05*i,'accepted')
                        else:
                            plt.text(lim1,0.05*i,'rejected')
                    plt.legend()                
            
            snr_filt = np.zeros(len(intervals)-1)
            if snr_filter_frequency_domain and snr_filter_time_domain:
                for i in range(len(snr_filt_fd)):
                    snr_filt[i] = snr_filt_fd[i] * snr_filt_td[i]
            elif snr_filter_frequency_domain:
                snr_filt = snr_filt_fd
            elif snr_filter_time_domain:
                snr_filt = snr_filt_td
     
    """ getting zero crossings and extracting phase velocity """    

    # crossings,phase_vel = noise.get_smooth_pv(
    #     ccdict['freq'],ccdict['spectrum'],
    #     dist,ref_curve,freqmin=min_freq,freqmax=max_freq,
    #     min_vel=min_vel, max_vel=max_vel,filt_width=3,filt_height=1.2,
    #     pick_threshold=pick_thresh,horizontal_polarization=horizontal_polarization,
    #     check_cross_quality=True,
    #     smooth_spectrum=True,plotting=plotresult)
                
    crossings,phase_vel = noise.get_smooth_pv(
        freqlist,speclist,
        dist,ref_curve,freqmin=min_freq,freqmax=max_freq,
        min_vel=min_vel, max_vel=max_vel,filt_width=3,filt_height=1.2,
        pick_threshold=pick_thresh,horizontal_polarization=horizontal_polarization,
        check_cross_quality=True,
        smooth_spectrum=not(velocity_filter),plotting=plotresult)
          
    pause
    if len(phase_vel)==0:
        accepted=False
        fails += 1
    else:
        accepted=True
        
    if snr_filter_frequency_domain or snr_filter_time_domain:
        freq_filtered = np.array([])
        phase_vel_filtered = np.array([])
        for i in range(len(intervals)-1):
            if snr_filt[i]:
                lim1 = intervals[i]
                lim2 = intervals[i+1]
                phase_vel_filtered = np.append(phase_vel_filtered,phase_vel[:,1][(lim1<=phase_vel[:,0])*(phase_vel[:,0]<lim2)])
                freq_filtered = np.append(freq_filtered,phase_vel[:,0][(lim1<=phase_vel[:,0])*(phase_vel[:,0]<lim2)])
        if statpair and subspec=='symmetric':
            plt.figure()
            plt.plot(1./phase_vel[:,0],phase_vel[:,1],'o',ms=6,label='without snr filter')
            plt.plot(1./freq_filtered,phase_vel_filtered,'o',ms=4,label='with snr filter')
            plt.legend()
        phase_vel = np.column_stack((freq_filtered,phase_vel_filtered))
                
    outfilename = os.path.basename(fname)
    outfilename = outfilename.replace("."+outfilename.split(".")[-1],".curve")
    if accepted and len(phase_vel)>0:
        np.savetxt(os.path.join(phase_vel_path,outfilename),phase_vel,
                   fmt="%.5f %.3f")
    elif len(phase_vel)>0:
        np.savetxt(os.path.join(phase_vel_path_lq,outfilename),phase_vel,
                   fmt="%.5f %.3f")
        
    counter += 1
    
    if statpair or plotresult or counter%500 == 0:
        plt.ioff()
        fig = plt.figure(figsize=(16,8))
        plt.suptitle(stat1+"_"+stat2+"   dist: %dkm   SNR: %.1f" %(dist,overall_snr))
        plt.subplot(221)
        plt.plot(timeax,np.fft.fftshift(tcorr),'k',linewidth=0.8,
                 label='raw cc (full stack)')
        if velocity_filter:
            plt.plot(timeax,np.fft.fftshift(np.fft.irfft(np.real(np.fft.rfft(TCORR_full)))),'r--',linewidth=0.6,
                     label='cc filtered and symmetric')
        plt.legend(numpoints=1)
        plt.xlim(dist/min_vel*-1.3,dist/min_vel*1.3)
        plt.title("CC")
        plt.xlabel('time lag [s]')
        plt.subplot(222)
        plt.plot(ccdict['freq'],np.zeros(len(ccdict['freq'])),linewidth=0.8,
                 color='darkgrey')
        plt.plot(ccdict['freq'],ccdict['spectrum'].real,'k',linewidth=1.2,
                 label='cc raw (full stack)')
        if velocity_filter:
            plt.plot(ccdict['freq'],np.real(np.fft.rfft(TCORR_full)),'r--',
                     linewidth=1,label='cc filtered')
        plt.plot(crossings[:,0],np.zeros(len(crossings)),'o',markeredgecolor='darkgrey',
                 linewidth=0.1,markerfacecolor='None',label='zero crossings',ms=3)
        plt.legend(numpoints=1)
        plt.title("real CC spectrum")
        plt.xlabel('frequency [Hz]')
        plt.xlim(min_freq,max_freq)
        ax = plt.subplot(223)
        if len(phase_vel)>0:
            plt.plot(1./phase_vel[:,0],phase_vel[:,1],'ro',ms=4,label='picked phase vels')
        plt.plot(1./crossings[:,0],crossings[:,1],'ko',ms=1,label='crossings')
        plt.plot(1./ref_curve[:,0],ref_curve[:,1],label='reference curve')
        plt.legend(numpoints=1)
        plt.xlim(1./np.max(crossings[:,0]),1./min_freq)
        plt.gca().set_xscale('log')
        plt.xlabel('period [s]')
        plt.ylabel('velocity [km/s]')
        plt.subplot(224)
        if len(phase_vel)>0:
            plt.plot(1./phase_vel[:,0],phase_vel[:,1],'k-',
                     linewidth=2,label='final phase velocity')
        linewidth=1
        # for ic in range(len(phasevel_curves)):
        #     if ic == 0:
        #         label='symmetric'
        #         color='red'
        #     elif ic == 1:
        #         label='from pos. lag time'
        #         color='green'
        #     elif ic == 2:
        #         label='from neg. lag time'
        #         color='cyan'
        #     elif ic == 3:
        #         label='from monthly correlation'
        #         color='black'
        #         linewidth=0.5
        #     else:
        #         label='_'
        #     plt.plot(1./ref_curve[:,0],phasevel_curves[ic],'-',color=color,
        #              linewidth=linewidth,label=label)
        plt.legend()
        plt.xlabel('period [s]')
        plt.ylabel('velocity [km/s]')
        if snr_filter_frequency_domain or snr_filter_time_domain:
            y0=np.min(phase_vel[:,1])-0.2
            for i in range(len(intervals)-1):
                lim1 = intervals[i]
                lim2 = intervals[i+1]
                if lim1==0:
                    lim1 = 1./300.
                if snr_filt[i]:
                    plt.plot([1./lim2,1./lim1],[y0,y0],color='green',lw=3)
                else:
                    plt.plot([1./lim2,1./lim1],[y0,y0],color='red',lw=3)
            plt.plot([],[],lw=3,color='red',label='snr filter: discarded')
            plt.plot([],[],lw=3,color='green',label='snr filter: accepted')
            plt.legend(numpoints=1)
        plt.gca().set_xscale('log')
        plt.xlim(1./np.max(crossings[:,0]),100)
        plt.ylim(min_vel,max_vel)
        if accepted:
            if statpair or plotresult:
                print("accepted")
            plt.savefig(figure_path+stat1+"_"+stat2+".png",bbox_inches='tight')
        else:
            print("quality check failed")
            plt.savefig(figure_path_failed+stat1+"_"+stat2+".png",bbox_inches='tight')
        if statpair:
            plt.show()
        plt.close(fig)


#%%       
#        plt.figure(figsize=(16,8))
#        plt.suptitle(stat1+"_"+stat2+"   dist: %dkm   SNR: %.1f" %(dist,overall_snr))
#        plt.subplot(221)
#        plt.plot((np.arange(len(spec_symm)*2-2)-len(spec_symm)+1)*dt,np.fft.fftshift(np.fft.irfft(spec_symm)),label='raw cc')
#        if velocity_filter:
#            plt.plot((np.arange(len(spec_symm)*2-2)-len(spec_symm)+1)*dt,np.fft.fftshift(tcorr_symm),'--',label='cc filtered and symmetric')
#        plt.legend(numpoints=1)
#        plt.xlim(dist/min_vel*-2.,dist/min_vel*2)
#        plt.title("CC")
#        plt.xlabel('time lag [s]')
#        plt.subplot(222)
#        plt.plot(freqax,np.real(spec_symm),label='cc raw')
#        if velocity_filter:
#            plt.plot(freqax,np.real(corr_symm),'--',label='cc filtered')
#        plt.legend(numpoints=1)
#        plt.title("real CC spectrum")
#        plt.xlabel('frequency [Hz]')
#        plt.subplot(223)
#        plt.plot(1./phase_vel_final[:,0],phase_vel_final[:,1],'ro',ms=4,label='picked phase vels')
#        plt.plot(1./crossings_symm[:,0],crossings_symm[:,1],'o',ms=2,label='crossings')
#        plt.plot(1./ref_curve[:,0],ref_curve[:,1],label='reference curve')
#        plt.legend(numpoints=1)
#        plt.xlim(min_period,50)
#        plt.xlabel('period [s]')
#        plt.ylabel('velocity [km/s]')
#        plt.subplot(224)
#        plt.plot(1./phase_vel_final[:,0],phase_vel_final[:,1])
#        plt.xlabel('period [s]')
#        plt.ylabel('velocity [km/s]')
#        if snr_filter_frequency_domain or snr_filter_time_domain:
#            y0=np.min(phase_vel_final[:,1])-0.2
#            for i in range(len(intervals)-1):
#                lim1 = intervals[i]
#                lim2 = intervals[i+1]
#                if lim1==0:
#                    lim1 = 1./300.
#                if snr_filt[i]:
#                    plt.plot([1./lim2,1./lim1],[y0,y0],color='green',lw=3)
#                else:
#                    plt.plot([1./lim2,1./lim1],[y0,y0],color='red',lw=3)
#            plt.plot([],[],lw=3,color='red',label='snr filter: discarded')
#            plt.plot([],[],lw=3,color='green',label='snr filter: accepted')
#            plt.legend(numpoints=1)
#        plt.xlim(3,100)
#        plt.show()
       
#np.save("fw_bw_diffcollection_rank%d_zz.npy" %mpi_rank,np.array(quality_check_collection))
