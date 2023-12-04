#import obspy
from obspy.clients.fdsn.mass_downloader import RectangularDomain, \
    Restrictions, MassDownloader
from obspy.clients.fdsn import Client
#from obspy.clients.fdsn import RoutingClient
from obspy import UTCDateTime
from obspy import read_inventory
import os
import sys
import errno
import logging
logger = logging.getLogger("obspy.clients.fdsn.mass_downloader")
logger.setLevel(logging.DEBUG)


def get_mseed_storage(network, station, location, channel, starttime, endtime):
    # Returning True means that neither the data nor the StationXML file
    # will be downloaded.
    filepath = os.path.join( "/media/HD3/database_alparray/%s/%s/%s/%s.D/%s.%s.%s.%s.D.%s.%03d" 
                             % (starttime.year,network,station,channel,network, station, 
                                location, channel,starttime.year,starttime.julday))

    # check also if this file has already been downloaded for another channel (HH, BH, CH)
    filepaths_alternative = []
    for cha_code in ["HH","BH","CH"]:
        filepaths_alternative.append(os.path.join( "/media/HD3/database_alparray/%s/%s/%s/%s.D/%s.%s.%s.%s.D.%s.%03d" 
                             % (starttime.year,network,station,cha_code+channel[2:],network, station, 
                                location, cha_code+channel[2:],starttime.year,starttime.julday)))

    skip_processed = False
    if skip_processed:
    # if True is returned, the channel is completely ignored. If a filename is
    # returned, it checks whether the file already exists and if yes, only the
    # station metadata is downloaded.
        if (os.path.isfile(filepath) or 
            os.path.isfile(filepaths_alternative[0]) or 
            os.path.isfile(filepaths_alternative[1]) or 
            os.path.isfile(filepaths_alternative[2]) ):
            #print("file already processed",filepath)
            return True
    
    #if not starttime.julday == 358:
    #    return True
    
    filepath = filepath.replace("database_alparray","SwathD")
    if os.path.isfile(filepath):
        return True
    
    return filepath
    
    #if is_in_db(network, station, location, channel, starttime, endtime):
    #    return True
    ## If a string is returned the file will be saved in that location.
    ##return os.path.join( "Albania_M6_MD/mseed/%s.%s.%s.%s.D" % (network, station, location, channel)
    #return 


#origin_time = UTCDateTime(2019, 11, 26, 2, 50, 00)

# Circular domain around the epicenter. This will download all data between
# 70 and 90 degrees distance from the epicenter. This module also offers
# rectangular and global domains. More complex domains can be defined by
# inheriting from the Domain class.

domain = RectangularDomain(minlatitude=44.5, maxlatitude=49,  
                           minlongitude=8, maxlongitude=18)

#inv_name= "%s" % sys.argv[1]
#inv=read_inventory(inv_name)
#print (inv)

restrictions = Restrictions(
    # Get data from 5 minutes before the event to one hour after the
    # event. This defines the temporal bounds of the waveform data.
    starttime=UTCDateTime(2017, 8, 1), # (2017, 8, 1) Swath D stations were in the field from Aug. 2017
    endtime=UTCDateTime(2020, 1, 1), #2020 1 1
    chunklength_in_sec=86400,
    # You might not want to deal with gaps in the data. If this setting is
    # True, any trace with a gap/overlap will be discarded.
    reject_channels_with_gaps=False,
    # And you might only want waveforms that have data for at least 95 % of
    # the requested time span. Any trace that is shorter than 95 % of the
    # desired total duration will be discarded.
    minimum_length=0.0,
    # if sanitize = True, stations with no station information are deleted
    sanitize = False,
    # No two stations should be closer than 10 km to each other. This is
    # useful to for example filter out stations that are part of different
    # networks but at the same physical station. Settings this option to
    # zero or None will disable that filtering.
    #network="ZS",
    #exclude_networks=["Z3"],
    #station="D104",
    minimum_interstation_distance_in_m=0,
    # Only HH or BH channels. If a station has HH channels, those will be
    # downloaded, otherwise the BH. Nothing will be downloaded if it has
    # neither. You can add more/less patterns if you like.
    channel_priorities=["BHZ", "HHZ"],#"HH[ZEN123]","CH[ZEN123]"],
    # Location codes are arbitrary and there is no rule as to which
    # location is best. Same logic as for the previous setting.
    location_priorities=["", "00", "01","02"],)
    #limit_stations_to_inventory=inv)


# No specified providers will result in all known ones being queried.
#clientroutingEIDA = RoutingClient('eida-routing', credentials={'EIDA_TOKEN': '/Users/irene/.eidatoken'})
ODC = Client('ODC', eida_token='/home/emanuel/.eidatoken',timeout=180)
#ODC = Client('ODC', timeout=180)
BGR = Client('BGR', timeout=180)
LMU = Client('LMU', eida_token='/home/emanuel/.eidatoken',timeout=180)
ETH = Client('ETH', eida_token='/home/emanuel/.eidatoken',timeout=180)
#NEIP = Client('NEIP')
GFZ = Client('GFZ', eida_token='/home/emanuel/.eidatoken',timeout=180)
NOA = Client('NOA',timeout=180)
#INGV = Client('INGV')
INGV = Client('INGV', eida_token='/home/emanuel/.eidatoken',timeout=180)
#KOERI = Client('KOERI',timeout=180)
#RESIF = Client('RESIF', user='494286_ingv')
RESIF = Client('RESIF',eida_token='/home/emanuel/.eidatoken',timeout=180)
IRIS = Client('IRIS',timeout=180)

mdl = MassDownloader(providers=[ BGR, LMU, ETH, GFZ, INGV, NOA, RESIF, ODC, IRIS ])
#mdl = MassDownloader(providers=[ ODC, IRIS ])

# The data will be downloaded to the ``./waveforms/`` and ``./stations/``
# folders with automatically chosen file names.
mseed_storage = get_mseed_storage 
stationxml_storage = "/media/HD3/SwathD/station_inventory/"
mdl.download(domain, restrictions, mseed_storage=mseed_storage,
             stationxml_storage=stationxml_storage, download_chunk_size_in_mb=50,
             threads_per_client=2)
