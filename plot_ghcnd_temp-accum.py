#!/usr/bin/python

###############################################################################
# plot_ghcnd_temp-accum.py
#
# Written by Jared Rennie (@jjrennie) 
#
# **DESCRIPTION**
# This program will take one station from GHCN-D and plot accumulated days 
# That are greater than some temperature threshold.
#
# **COPYRIGHT**
# THIS SOFTWARE AND ITS DOCUMENTATION ARE CONSIDERED TO BE IN THE PUBLIC
# DOMAIN AND THUS ARE AVAILABLE FOR UNRESTRICTED PUBLIC USE.  THEY ARE
# FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES GOVERNMENT, ITS
# INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND AGENTS MAKE NO WARRANTY,
# EXPRESS OR IMPLIED, AS TO THE USEFULNESS OF THE SOFTWARE AND
# DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME NO RESPONSIBILITY (1) FOR
# THE USE OF THE SOFTWARE AND DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL
# SUPPORT TO USERS.
###############################################################################

# Import Modules
import numpy as np
import numpy.ma as ma
import sys 
import time
import datetime
import calendar
import re
import pylab
from calendar import monthrange
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from datetime import datetime, date
import pandas as pd
from ftplib import FTP

# Declare User Directories
main_directory="."

# Declare Other Variables
begin_year=1700
num_elements=1 
tmax=0
num_days=366
normals_start=1991
normals_end=2020

# Begin Program
start=time.time()

# Read in Arguments 
if len(sys.argv) < 3:
    sys.exit("USAGE: <GHCND ID> <TEMP_THRESH>")  
station = sys.argv[1]
temp_thresh=float(sys.argv[2])
print(station,str(temp_thresh))

end_year=datetime.now().year 
num_years=(end_year-begin_year) + 1

#################################################
# Read in GHCND-D Meta 
print("\nGETTING METADATA FOR STATION: ",station)

ftp = FTP('ftp.ncei.noaa.gov')
ftp.login()
ftp.cwd('pub/data/ghcn/daily')
ftp.retrbinary('RETR '+'ghcnd-stations.txt', open('ghcnd-stations.txt', 'wb').write)
ftp.quit()

ghcnd_stnfile=main_directory+'/ghcnd-stations.txt'
ghcnd_stations= np.genfromtxt(ghcnd_stnfile,delimiter=(11,9,10,7,4,30),dtype=str)

ghcnd_meta = ghcnd_stations[ghcnd_stations[:,0] == station]

ghcnd_id=ghcnd_meta[0][0]
ghcnd_lat=float(ghcnd_meta[0][1])
ghcnd_lon=float(ghcnd_meta[0][2])
ghcnd_alt=float(ghcnd_meta[0][3])
ghcnd_state=ghcnd_meta[0][4]
ghcnd_name=ghcnd_meta[0][5]
ghcnd_name = ghcnd_name.strip()
ghcnd_name = re.sub(' +',' ',ghcnd_name)
ghcnd_name = ghcnd_name.replace(" ","_")

#################################################
# Read in GHCN-D Data (Original, QC'd data removed)
print("\nGETTING DATA FOR STATION: ",station)

ftp = FTP('ftp.ncei.noaa.gov')
ftp.login()
ftp.cwd('pub/data/ghcn/daily/all')
ftp.retrbinary('RETR '+ghcnd_id+'.dly', open(ghcnd_id+'.dly', 'wb').write)
ftp.quit()

infile = main_directory+'/'+ghcnd_id+".dly"
ghcnd_value = np.zeros((num_years,12,31,num_elements),dtype='f')
tmax_nonmiss = np.zeros((num_years,num_elements),dtype='f')

file_handle = open(infile, 'r')
contents = file_handle.readlines()
file_handle.close() 

valid_end=-9999
valid_begin=9999
num_normal=0
for counter in range(len(contents)): 
    element = contents[counter][17:21]
    if element == "TMAX":
        element_counter=tmax
        year = int(contents[counter][11:15])
        if year <= end_year:
            year_counter = year-begin_year
            valid_begin=min(valid_begin,year)
            valid_end=max(valid_end,year)
            month = int(contents[counter][15:17])
            month_counter = month-1
            char=21
            for day_counter in range(0,31):
                if contents[counter][char:char+5] != "-9999" and contents[counter][char+6:char+7] == " ": # Remove QC
                    ghcnd_value[year_counter][month_counter][day_counter][element_counter] = float(contents[counter][char:char+5]) / 10.0
                    tmax_nonmiss[year_counter][element_counter]+=1
                    if year >= normals_start and year <= normals_end:
                        num_normal+=1
                    last_day=day_counter+1
                char = char + 8

# Get day of year for last day with valid data
last_date=datetime(year, month, last_day)
last_day=datetime(year, month, last_day).timetuple().tm_yday
valid_years=(valid_end-valid_begin) + 1

# Convert from C to F
ghcnd_value=(ghcnd_value*1.8) + 32

#################################################
# Create Accumulations
raw_tmax = np.zeros((num_years,num_days),dtype='f')
tmax_accum = np.zeros((num_years,num_days),dtype='f')
total_accum = np.zeros((num_years),dtype='f')

for year_counter in range(0,num_years):
    day_of_year=0
    day_before=0
    for month_counter in range(0,12):
        for day_counter in range(0,31):
            try:
                # Check if date is valid
                datetime(year=year_counter+begin_year,month=month_counter+1,day=day_counter+1)
                raw_tmax[year_counter][day_of_year] = ghcnd_value[year_counter,month_counter,day_counter,tmax]
                if ghcnd_value[year_counter,month_counter,day_counter,tmax] >= temp_thresh:
                    tmax_accum[year_counter][day_of_year] = day_before + 1
                else:
                    tmax_accum[year_counter][day_of_year] = day_before
                total_accum[year_counter]=tmax_accum[year_counter][day_of_year]
                day_before=tmax_accum[year_counter][day_of_year]
                day_of_year=day_of_year+1
            except:
                pass   

# Get 30 Year Average (Straight From Data)
average_tmax = np.zeros((num_days),dtype='f')-(9999.0)
for day_counter in range(0,num_days):
    average_tmax[day_counter]=ma.average(tmax_accum[(normals_start-begin_year):(normals_end-begin_year),day_counter])
average_tmax=np.ceil(average_tmax)

#################################################
# PLOT

# Mask Zero Data before plotting
tmax_accum = ma.masked_values(tmax_accum, 0.)
total_accum = ma.masked_values(total_accum, 0.)

#Get Some Stats Needed For Plotting
x_axis=range(num_days)
x_axis_end=range(last_day)

# Current Year
current_loc = num_years-1
current_tmax = "%3i" % total_accum[current_loc]
current_year = current_loc + begin_year
current_data=tmax_accum[current_loc,0:last_day]
current_last=tmax_accum[current_loc,last_day]

max_tmax = "%3i" % np.max(total_accum)
max_loc = np.argmax(total_accum)
max_year = max_loc+begin_year

min_tmax = "%3i" % np.min(total_accum[np.where(total_accum != 0)])
min_loc = np.nanargmin(total_accum)
min_year = min_loc+begin_year

#Need complete year for Min Value
counter=0
not_full_year=True
sorted=np.argsort(total_accum)
while not_full_year:
    if tmax_nonmiss[min_loc][0] < 300:
        counter+=1
        min_loc=sorted[counter]
    else:
        min_year = min_loc+begin_year
        min_tmax = "%3i" % np.ma.max(tmax_accum[min_loc,:])
        not_full_year=False

# Average Year
avg_tmax = "%3i" % average_tmax[365]
if average_tmax[-1] < average_tmax[-2]:
    avg_tmax = "%3i" % average_tmax[-2]

# Create Figure
fig, ax_in = plt.subplots(figsize=(15, 8), edgecolor='white', facecolor='white', dpi=300)
plt.style.use("dark_background")

# Add grid lines
plt.grid(color='white', linestyle='--', linewidth=0.5, alpha=0.3)
ax_in.set_facecolor('#808080')

# Plot Accumulated Val (Sort by end of year accumulation and plot by range of color)
order=np.argsort(tmax_accum[:,364])
color_pos=np.linspace(0.5,1,valid_years)
order_counter=0
color_counter=0
for year_counter in range(0,num_years):
    pos=order[order_counter]
  
    if pos != (num_years-1) and not np.ma.is_masked(np.sum(tmax_accum[pos,:])):
        plt.plot(x_axis, tmax_accum[pos,:], linewidth=0.5, color=colors.rgb2hex(pylab.cm.YlOrRd(color_pos[color_counter])[0:3]))  
        color_counter=color_counter+1
    order_counter=order_counter+1

# Overlay Record Max Year
outHex=colors.rgb2hex(pylab.cm.YlOrRd(color_pos[len(color_pos)-1])[0:3]) 
if max_loc==current_loc:
    plt.plot(x_axis_end, tmax_accum[max_loc,0:last_day], color=outHex, linewidth=3, label='Max ('+str(max_year)+': '+str(max_tmax)+')')  
else:
    plt.plot(x_axis, tmax_accum[max_loc,:], color=outHex, linewidth=3, label='Max ('+str(max_year)+': '+str(max_tmax)+')')  
  
# Overlay Record Min Year
outHex=colors.rgb2hex(pylab.cm.YlOrRd(color_pos[0])[0:3]) 
if min_loc==current_loc:
    plt.plot(x_axis_end, tmax_accum[min_loc,0:last_day], color=outHex, linewidth=3, label='Min ('+str(min_year)+': '+str(min_tmax)+')')
else:
    plt.plot(x_axis, tmax_accum[min_loc,:], color=outHex, linewidth=3, label='Min ('+str(min_year)+': '+str(min_tmax)+')')  

# Overlay Average 
plt.plot(x_axis[0:365], average_tmax[0:365], color='#e6b800', linewidth=3, markeredgecolor='white', label='Avg ('+str(normals_start)+'-'+str(normals_end)+': '+str(avg_tmax)+')') 

# Overlay Current Year
plt.plot(x_axis_end, current_data, color='black', markeredgecolor='white', linewidth=3, label='Current ('+str(current_year)+': '+str(current_tmax)+')')    
plt.plot(x_axis_end[last_day-1],current_last, marker='o', color='black', markersize=10)

# Plot Legend
plt.legend(bbox_to_anchor=(0., -.102, 1., -1.02), loc=3, ncol=4, mode="expand", borderaxespad=0., fontsize=11, facecolor='#808080')

# Plot X/Y Limits
ymin=0
ymax=int(5 * round(float((np.max(tmax_accum) + 10))/5))
plt.ylim(ymin,ymax)
plt.xlim(-5, num_days) 

# Plot Y-Axis Label
plt.yticks(range(ymin, ymax, 10), [r'{}'.format(x) for x in range(ymin, ymax, 10)], fontsize=10, color='white')
plt.ylabel(r'Number of Accumulated Days', fontsize=12, color='white')

# Plot X-Axis Label
month_pos=[0,31,60,91,121,152,182,213,244,274,305,335]
month_names=["Jan 1","Feb 1","Mar 1","Apr 1","May 1","Jun 1","Jul 1","Aug 1","Sep 1","Oct 1","Nov 1","Dec 1"]
plt.xticks(month_pos, month_names, fontsize=10, color='white')

# Plot 2nd Y Axis Labels
ax_mm = ax_in.twinx()
y1, y2 = ax_in.get_ylim()
ax_mm.set_ylim(int(y1), int(y2))
ax_mm.figure.canvas.draw()
ax_in.callbacks.connect("ylim_changed", ax_mm)
ax_mm.set_ylabel(r' Number of Accumulated Days', fontsize=12, rotation=270, labelpad=20)

# Get Day Percentage
valid_days=tmax_nonmiss[current_loc][0]
all_days = datetime.now().timetuple().tm_yday
day_percentage= (float(float(valid_days) / float(all_days)) * 100.)

# Plot Title/Subtitle
plt.suptitle(str(current_year)+' Number Days >= '+str(int(temp_thresh))+'F For '+ghcnd_name.title()+', '+ghcnd_state.strip(), fontsize=20)
plt.title('GHCN-D ID= '+station+' | LAT= '+str(ghcnd_lat)+' | LON= '+str(ghcnd_lon)+' | ELEV= '+str(int(ghcnd_alt*3.2808399))+'\' | Period= '+str(valid_begin)+'-'+str(valid_end)+' |  YTD Report= '+str(int(day_percentage))+'%', fontsize=12)
plt.annotate('Generated by Jared Rennie (@jjrennie) on: '+str(datetime.now())[0:16]+' Eastern. Data up to '+(last_date.strftime('%Y%m%d')),xy=(0.995, 0.01), xycoords='axes fraction', fontsize=7,horizontalalignment='right', verticalalignment='bottom')

# Save Figure
plt.savefig(main_directory+'/'+station+'_'+str(int(temp_thresh))+'F.png', dpi=300,bbox_inches='tight')
plt.clf()

####################
# DONE
####################
print("DONE!")
end=time.time()
print ("Runtime: %8.1f seconds." % (end-start)) 
sys.exit()
