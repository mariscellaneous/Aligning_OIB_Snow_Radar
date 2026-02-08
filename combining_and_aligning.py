import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import netCDF4 as nc
import glob
from scipy.io import loadmat
from matplotlib.image import imread
import matplotlib.pyplot as plt
import os
from scipy import interpolate
import datetime
import scipy 
import pickle
import cartopy.crs as ccrs
import scipy.spatial as spatial

def lining_up_lines(finished_file):
    plt.ioff()
    first_file = finished_file[0]
    
    # first_year = first_file[37:41]
    # first_month = first_file[41:43]
    # first_day = first_file[43:45]
    
    # first_date = datetime.datetime(int(first_year),int(first_month),int(first_day))
    path =  first_file[:76]+'/close_files'
    pathext= path+'/'+first_file[-27:-4]
    
    
    
    # isExist = os.path.isfile(first_file[:46]+'close_files_pic/'+first_file[-27:-4]+'_NearbyFiles_with200_3_new.png')
    # if isExist==True:
    #     return'
    second_file_temp=finished_file[-1]
    isExist2 = os.path.isfile(path+'/'+first_file[-27:-4]+'_'+second_file_temp[-27:-4]+'matching_file_FINAL.png')
    # if isExist2==True:
    #     return
    
    first_ds = loadmat(first_file)
    
    first_lat = first_ds['Latitude'][0]
    first_lon = first_ds['Longitude'][0]
    
    
    
    first_time_orig = np.squeeze(np.arange(0,len(first_ds['Data'][200:,0]))*first_ds['del_t']*10**-6)
    
    
    first_time = np.squeeze(np.arange(0,len(first_ds['Data'][200:,0]))*first_ds['del_t']*10**-6)
    first_data = first_ds['Data'][200:,:]

    
    next_file_ind = 1
    second_dats = []
    second_ts = []
    second_ds = []
    
    for next_file_ind,second_file in enumerate(finished_file[1:]):
        path =  first_file[:46+34-4]+'/close_files'
        second_ds_mat = loadmat(second_file[:-4])
        second_lat = second_ds_mat['Latitude'][0]
        second_lon = second_ds_mat['Longitude'][0]
        second_time = np.squeeze(np.arange(0,len(second_ds_mat['Data'][200:,0]))*second_ds_mat['del_t']*10**-6)
    
    
        second_data = np.zeros((len(second_ds_mat['Data'][0,:]),len(first_time))).T*np.nan
    
        for along_data_ind in np.arange(0, second_data.shape[1]):
            y = second_ds_mat['Data'][200:,along_data_ind]
            f = interpolate.interp1d(np.squeeze(second_time),y,fill_value="extrapolate")
            second_data[:,along_data_ind]=f(first_time)
        
        maxarraysize=np.max([first_data.shape[0],second_data.shape[0]])
        d,t = np.meshgrid(np.arange(0,first_data.shape[1])/10.,first_time)
        d2,t2 = np.meshgrid(np.arange(0,first_data.shape[1])/10.,first_time)
        
    
        
    
    
        allocate_first = np.zeros((maxarraysize,first_data.shape[1])).T*np.nan
        allocate_second = np.zeros((maxarraysize,second_data.shape[1])).T*np.nan
    
        allocate_first[:first_data.shape[1],:first_data.shape[0]] = first_data.T.copy()
        
        allocate_second[:second_data.shape[1],:second_data.shape[0]] = second_data.T.copy()
    
    
              
    
    
        alldata = np.hstack((allocate_first.T,allocate_second.T)).T
        
        mean200_pri = np.nanmean(first_ds['Data'][:200,:])
        mean200_sec = np.nanmean(second_ds_mat['Data'][:200,:])
    
        stdev200_pri = np.mean(np.nanstd(first_ds['Data'][:200,:],1))
        stdev200_sec = np.mean(np.nanstd(second_ds_mat['Data'][:200,:],1))
    
        # alldata = np.hstack((first_data,second_data))
    
        new_col = alldata.sum(1)[...,None] # None keeps (n, 1) shape
        new_col[:,0] = np.nan
        alldata = np.append(alldata, new_col, 1)
    
    
    
        zee_point = projection.transform_points(ccrs.PlateCarree(),np.asarray(first_lon),np.asarray(first_lat))
        zee_point_2 = projection.transform_points(ccrs.PlateCarree(),np.asarray(second_lon),np.asarray(second_lat))
    
        xvals_first = zee_point[:,0]
        yvals_first = zee_point[:,1]
        xvals_second = zee_point_2[:,0]
        yvals_second = zee_point_2[:,1]
        spacing = np.sqrt((xvals_first[2]-xvals_first[1])**2 + (yvals_first[2]-yvals_first[1])**2)

        points = np.vstack((xvals_second,yvals_second)).T
    
        point_tree2 = spatial.KDTree(points)
    
        cats = np.zeros((len(xvals_first),20))-1
        distances=[]
        for num in np.arange(0,len(xvals_first)):
            cat2 = point_tree2.query([xvals_first[num],yvals_first[num]],k=50,)
            cat = cat2[1]
            max_distance = 100
            cat[cat2[0]>max_distance] = -1
            distances.append(cat2[0][0])
            cats[num] = int(cat[0]) 
    
    
            
                
    
        del cat2, point_tree2
    
    
        if next_file_ind == 0:
            first_t = t[200:,:len(first_lat)]-t[200,:len(first_lat)]
            first_d = d[200:,:len(first_lat)]-d[200,:len(first_lat)]
            indices = cats.astype(int)
            zee_point1 = projection.transform_points(ccrs.PlateCarree(),np.asarray(first_lon),np.asarray(first_lat))

        
    
        zee_point2 = projection.transform_points(ccrs.PlateCarree(),np.asarray(second_lon),np.asarray(second_lat))
        
        indices = cats[:,1].astype(int)
    
        new_col = allocate_second.T.sum(1)[...,None] # None keeps (n, 1) shape
        new_col[:,0] = np.nan
        allocate_second_filled = np.append(allocate_second.T, new_col, 1).T
        
    
        new_data = np.zeros(allocate_first.shape)
    
        for ind,indice in enumerate(indices):
            new_data[ind,:] = allocate_second_filled[indice,:]
            
            if indice == -1:
                continue
    
        second_dat = new_data
    
    
        if np.sum(~np.isnan(new_data))==0:
            continue
    
    
    
    
        fig_extra, (ax_extra1,ax_extra2) = plt.subplots(nrows=2, ncols=1, figsize=(6,18))
        
        ax_extra1.set_title(first_file[-27:-4])
        ax_extra2.set_title(second_file[-27:-4])
        indices = cats.astype(int)
    
    
    
        
        first_dat = first_data
    
        ax_extra1.pcolormesh(d,t,allocate_first.T,cmap='Blues',shading='nearest',)
        first_t = t[:,:len(first_lat)]
    
        ax_extra1.set_ylim(0.2*10**-6,0)
        ax_extra2.set_ylim(0.2*10**-6,0)
        ax_extra2.set_xlabel('Distance [km]')
        ax_extra2.set_ylabel('twtt [s]')
        ax_extra2.pcolormesh(d2,t2,new_data.T,cmap='Blues')
        
        isExist2 = os.path.exists(path)
        
        if isExist2 == False:
            os.makedirs(path)
            
        fig_extra.savefig(path+'/'+first_file[-27:-4]+'_'+second_file[-27:-4]+'_initial.png')
    
    
        
        scipy.io.savemat(path+'/'+first_file[-27:-4]+'_'+second_file[-27:-4]+'initial.mat',{'Primary_Data':allocate_first.T,\
                                                                        'Primary_Depths':first_d,\
                                                                        'Primary_Lat':first_lat,'Primary_Lon':first_lon,\
                                                                        'Secondary_Data':new_data.T,\
                                                                        'Mean_200_primary':mean200_pri,'Mean_200_sec':mean200_sec,\
                                                                        'StDev_200_primary':stdev200_pri,'StDev_200_sec':stdev200_sec,\
                                                                        'del_t2':first_ds['del_t'],'del_t1':first_ds['del_t'],'distances':distances})    
       
    
        plt.tight_layout()
        
        isExist = os.path.exists(path)
    
        if isExist == False:
            os.makedirs(path)
    
    
    
    
    
    
    del first_t, first_d, second_dat, second_data, first_data, cats, distances, alldata, first_ds, second_ds_mat, allocate_first, allocate_second


def combining_files(finished_file):
    plt.ioff()
    first_file = finished_file[0]

    path =  first_file[:76]+'/close_files'
    pathext= path+'/'+first_file[-27:-4]
    
    first_ds = loadmat(first_file)
    
    first_lat = first_ds['Latitude'][0]
    first_lon = first_ds['Longitude'][0]


    
    excess_files = np.asarray(finished_file[1:])

    younger_files = np.asarray(excess_files)
    
    end_of_files = []
    begin_of_files = []
    for younger_file in younger_files:
        end_of_files.append(int(younger_file[-7:-4]))
        begin_of_files.append(younger_file[:-8])
    
    end_of_files=np.asarray(end_of_files)
    begin_of_files=np.asarray(begin_of_files)
    list_of_files = []
    for dates in excess_files:
        theslice =  np.abs(end_of_files - int(dates[-7:-4]))<2
        list_of_files.append(list(excess_files[(begin_of_files == dates[:-8]) & (theslice)]))
    


    sets_of_files = set(frozenset(i) for i in list_of_files)
    
    list_of_files = list(sets_of_files)

    new_list_of_files=[]
    for index,item in enumerate(sets_of_files):
        mysum=0
        for itemz in sets_of_files:
            if set(item).issubset(itemz):
                if item != itemz:
                    mysum+=1 
        if mysum == 0:
            new_list_of_files.append(item)
    
    lengths=[]

    for filelist in new_list_of_files:
        lengths.append(len(filelist))
    maxlength = max(lengths)
    list_of_files_padded = []
    str_len = len(first_file)
    for filelist in new_list_of_files:
        the_add = maxlength-len(filelist)
        tempfileset = list(filelist)
        if the_add > 0:
            for ind in np.arange(0,the_add):
                tempfileset.append("0"*str_len)
        list_of_files_padded.append(tempfileset)
    
    


    shorted_list_of_files = np.unique(list_of_files_padded,axis=0).tolist() # ADD THIS BACK IN!!

    if not isinstance(shorted_list_of_files[0], list):
        print('ooooh isinstance coming into play')
        shorted_list_of_files = [shorted_list_of_files]
    for file_item in shorted_list_of_files:
        temp_file_item = list(file_item)
        if "0"*112 in temp_file_item:
            temp_file_item.remove("0"*str_len)

        if len(temp_file_item) == 3:
            doesitexist_2sub1 = os.path.exists(path+'/'+first_file[-27:-4]+'_'+temp_file_item[0][-27:-4]+'initial.mat')
            doesitexist_2sub2 = os.path.exists(path+'/'+first_file[-27:-4]+'_'+temp_file_item[1][-27:-4]+'initial.mat')
            doesitexist_2sub3 = os.path.exists(path+'/'+first_file[-27:-4]+'_'+temp_file_item[2][-27:-4]+'initial.mat')

            if (doesitexist_2sub1==False) & (doesitexist_2sub2==False) & (doesitexist_2sub3==False):
                continue         

            flagging = np.asarray([doesitexist_2sub1,doesitexist_2sub2,doesitexist_2sub3])          
            temp_file_item_temp=np.asarray(temp_file_item)[flagging==True]
            temp_file_item = list(temp_file_item_temp)
            print(temp_file_item)

        if len(temp_file_item) == 1:

            doesitexist_1 = os.path.exists(path+'/'+first_file[-27:-4]+'_'+temp_file_item[0][-27:-4]+'initial.mat')
            if doesitexist_1 == False:
                continue                
        
        if len(temp_file_item) == 2:
            doesitexist_2sub1 = os.path.exists(path+'/'+first_file[-27:-4]+'_'+temp_file_item[0][-27:-4]+'initial.mat')
            doesitexist_2sub2 = os.path.exists(path+'/'+first_file[-27:-4]+'_'+temp_file_item[1][-27:-4]+'initial.mat')

            if (doesitexist_2sub1==False) & (doesitexist_2sub2==False):
                continue         
                
            flagging = np.asarray([doesitexist_2sub1,doesitexist_2sub2])     
                

            temp_file_item_temp=np.asarray(temp_file_item)[flagging==True]

            flagging=np.asarray(flagging)
            temp_file_item_temp=np.asarray(temp_file_item)[flagging==1]
            temp_file_item = list(temp_file_item_temp)
            
        if len(temp_file_item) == 3:
            print('ayoo, there are three files')
            doesitexist_2sub1 = os.path.exists(path+'/'+first_file[-27:-4]+'_'+temp_file_item[0][-27:-4]+'initial.mat')
            doesitexist_2sub2 = os.path.exists(path+'/'+first_file[-27:-4]+'_'+temp_file_item[1][-27:-4]+'initial.mat')
            doesitexist_2sub3 = os.path.exists(path+'/'+first_file[-27:-4]+'_'+temp_file_item[2][-27:-4]+'initial.mat')

            if (doesitexist_2sub1==False) & (doesitexist_2sub2==False) & (doesitexist_2sub3==False):
                continue         
                
            flagging = np.asarray([doesitexist_2sub1,doesitexist_2sub2,doesitexist_2sub3])          

            temp_file_item_temp=np.asarray(temp_file_item)[flagging==True]

            flagging=np.asarray(flagging)
            temp_file_item_temp=np.asarray(temp_file_item)[flagging==1]
            temp_file_item = list(temp_file_item_temp)

            if len(temp_file_item) > 2:
                temp_file_item = np.sort(temp_file_item)[:2]
                
        
            
        if len(temp_file_item) == 1:
            mydata = scipy.io.loadmat(path+'/'+first_file[-27:-4]+'_'+temp_file_item[0][-27:-4]+'initial.mat')
            filename = path+'/'+first_file[-27:-4]+'_'+temp_file_item[0][-27:-4]+'_add_none_new_FINAL.mat'

            if np.sum(~np.isnan(mydata['Secondary_Data'])) == 0:
                print('One file, and not enough data')
                continue
            scipy.io.savemat(filename,{'Primary_Data':mydata['Primary_Data'],\
                    'Primary_Lat':mydata['Primary_Lat'],'Primary_Lon':mydata['Primary_Lon'],\
                    'Secondary_Data':mydata['Secondary_Data'],\
                    'Mean_200_primary':mydata['Mean_200_primary'],'Mean_200_sec':mydata['Mean_200_sec'],\
                    'StDev_200_primary':mydata['StDev_200_primary'],'StDev_200_sec':mydata['StDev_200_sec'],\
                    'del_t2':mydata['del_t2'],'del_t1':mydata['del_t1']})   

            fig_extra, (ax_extra1,ax_extra2) = plt.subplots(nrows=2, ncols=1, figsize=(6,18))
            
            ax_extra1.set_title(first_file[-27:-4])
            ax_extra2.set_title(temp_file_item[0][-27:-4])

            [m,n]=np.meshgrid(np.arange(0,mydata['Primary_Data'].shape[1])*0.1,np.arange(0,mydata['Primary_Data'].shape[0]))
            ax_extra1.pcolor(m[::5,:],n[::5,:],mydata['Primary_Data'][::5,:],cmap='Blues')
            ax_extra2.pcolor(m[::5,:],n[::5,:],mydata['Secondary_Data'][::5,:],cmap='Blues')
            ax_extra1.invert_yaxis()
            ax_extra2.invert_yaxis()
            fig_extra.savefig(filename+'.png')
            print(filename+'.png')
            plt.close()
            del mydata
          
            
        if len(temp_file_item) == 2:
            mydata = scipy.io.loadmat(path+'/'+first_file[-27:-4]+'_'+temp_file_item[0][-27:-4]+'initial.mat')    
            mydata2 = scipy.io.loadmat(path+'/'+first_file[-27:-4]+'_'+temp_file_item[1][-27:-4]+'initial.mat')
            first_sec_data = mydata['Secondary_Data']
            first_sec_data[~np.isnan(mydata2['Secondary_Data'])] = np.nan
            sec_sec_data = mydata2['Secondary_Data']    
            total_sec_data=np.nansum(np.dstack((first_sec_data,sec_sec_data)),2)

            nanidx = np.isnan(first_sec_data)*np.isnan(sec_sec_data)

            
            total_sec_data[nanidx] = np.nan
            

            if np.sum(~np.isnan(total_sec_data)) == 0:
                print('Flights too far apart, no data here')
                continue
            
            filename = path+'/'+first_file[-27:-4]+'_'+temp_file_item[0][-27:-4]+'_add_'+temp_file_item[1][-27:-4]+'_FINAL.mat'
            scipy.io.savemat(filename,{'Primary_Data':mydata['Primary_Data'],\
                        'Primary_Lat':mydata['Primary_Lat'],'Primary_Lon':mydata['Primary_Lon'],\
                        'Secondary_Data':total_sec_data,\
                        'Mean_200_primary':(mydata['Mean_200_primary']+mydata2['Mean_200_primary'])/2,'Mean_200_sec':(mydata['Mean_200_sec']+mydata2['Mean_200_sec'])/2,\
                        'StDev_200_primary':(mydata['StDev_200_primary']+mydata2['StDev_200_primary'])/2,'StDev_200_sec':(mydata['StDev_200_sec']+mydata2['StDev_200_sec'])/2,\
                        'del_t2':mydata['del_t2'],'del_t1':mydata['del_t1']})     
            fig_extra, (ax_extra1,ax_extra2) = plt.subplots(nrows=2, ncols=1, figsize=(6,18))
            
            ax_extra1.set_title(first_file[-27:-4])
            ax_extra2.set_title(temp_file_item[0][-27:-4]+', '+temp_file_item[1][-27:-4])

            [m,n]=np.meshgrid(np.arange(0,mydata['Primary_Data'].shape[1])*0.1,np.arange(0,mydata['Primary_Data'].shape[0]))
            cax2=ax_extra1.pcolor(m[::5,:],n[::5,:],mydata['Primary_Data'][::5,:],cmap='Blues')
            cax=ax_extra2.pcolor(m[::5,:],n[::5,:],total_sec_data[::5,:],cmap='Blues')
            plt.colorbar(cax2,ax=ax_extra1)
            plt.colorbar(cax,ax=ax_extra2)
            ax_extra1.invert_yaxis()
            ax_extra2.invert_yaxis()
            fig_extra.savefig(filename+'.png')
            print(filename+'.png')
            plt.close()

            if len(temp_file_item)>2:
                print(temp_file_item)
                print(len(temp_file_item))
                print('OHHHH NOOO!!')
                
            del mydata, mydata2

            
    

    
