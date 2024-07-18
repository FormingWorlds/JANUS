import os
from osfclient.api import OSF

basic_list =[
        "/Dayspring/256",
        "/Frostflow/256",
        "/Legacy",
        "/Mallard",
        "/Oak",
        "/Reach"
        ]

def download_folder(storage, folder_name, local_path):
    ''''
    Download a specific folder in the OSF repository
    
    Inputs :
        - storage     : OSF storage name 
        - folder_name : folder name to be downloaded
        - local_path  : local repository where data are saved
    '''
    for file in storage.files:
        if file.path.startswith(folder_name):
            local_file_path = local_path + file.path
            #Create local directory if needed
            os.makedirs(os.path.dirname(local_file_path), exist_ok=True)
            #Download the file
            with open(local_file_path, 'wb') as local_file:
                file.write_to(local_file)
    return


def GetFWLData():
    fwl_data_dir = os.getenv('FWL_DATA')
    if os.environ.get("FWL_DATA") == None:
        raise Exception("The FWL_DATA environment variable where spectral data will be downloaded needs to be set up!")
    return fwl_data_dir

def DownloadStellarSpectra():
    #project ID of the stellar spectra on OSF 
    project_id = '8r2sw'

    # Link with OSF project repository
    osf = OSF()
    project = osf.project(project_id)
    storage = project.storage('osfstorage')

    # Folder
    data_dir = GetFWLData() + "/stellar_spectra"
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)

    # Get all named spectra
    download_folder(storage,"/Named",data_dir)


def DownloadSpectralFiles(fname="",nband=256):
    ''''
    Download spectral files data
    
    Inputs :
        - fname (optional) :    folder name, i.e. "/Dayspring"
                                if not provided download all the basic list  
        - nband (optional) :    number of band = 16, 48, 256, 4096 
                                (only relevant for Dayspring, Frostflow and Honeyside)
    '''

    #project ID of the spectral files on OSF 
    project_id = 'vehxg'

    #Create spectral file data repository if not existing
    data_dir = GetFWLData() + "/spectral_files"
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)

    #Link with OSF project repository
    osf = OSF()
    project = osf.project(project_id)
    storage = project.storage('osfstorage')

    #If no folder specified download all basic list
    if not fname:
        for folder in basic_list:
            if not os.path.exists(data_dir+folder):
                download_folder(storage,folder,data_dir)
    elif fname in ["/Dayspring","/Frostflow","/Honeyside"]:
        folder = fname + "/" + str(nband)
        if not os.path.exists(data_dir+folder):
            download_folder(storage,folder,data_dir)
    elif fname in ["/Kynesgrove","/Legacy","/Mallard","/Oak","/Reach","/stellar_spectra"]:
        folder = fname
        if not os.path.exists(data_dir+folder):
            download_folder(storage,folder,data_dir)
    else:
        print("Unrecognised folder name in DownloadSpectralFiles")

    return
