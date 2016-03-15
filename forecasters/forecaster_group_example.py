import time
import subprocess
import sys

def CheckHaltFile(haltfilename):
	halt = 0
	with open(haltfilename,'r') as haltfile:
		halt = int(haltfile.readline())
	return halt

np = sys.argv[1]

#Number of forecasters in the group
num_forecasters = 1

#Filenames
haltfilename = 'examples/fgroup/halt'
timesfilename = 'examples/fgroup/times'
ptimesfilename = 'examples/fgroup/output_'

#Create halt file
halt = 0
with open(haltfilename,'w') as haltfile:
	haltfile.write(str(halt))

#Get the first batch of init and final times
with open(timesfilename) as infile:
	old_times = []
	for line in infile:
		if line.strip():
			old_times.append([int(x) for x in line.split()])

#Check that enough initial timestamps are set in the times file
N = len(old_times)
if N != num_forecasters:
	print 'Error: expected',num_forecasters,'forecasters. Got data for ',N
	sys.exit(1)

while(halt == 0):

	#Call programs #######################
	idx = -1

	#0: Toplayer - IFC (ifc1c)
	idx += 1
	cmd = 'mpirun -np '+str(np)+' ./FORECASTER_MAPS_END examples/GlobalForecast262_ifc1c.gbl examples/fcast_file.fcst '+str(old_times[idx][0])+' '+str(old_times[idx][1])+' '+ptimesfilename+str(idx)+' '+str(old_times[idx][0]-3600)+' 0 0'
	print '\nRunning command',cmd
	sys.stdout.flush()
	time.sleep(1)
	subprocess.call(cmd,shell=True)

	#Create new times files
	new_times = []
	outfile = open(timesfilename,'w')
	for i in range(N):
		with open(ptimesfilename+str(i),'r') as infile:
			for line in infile:
				towrite = ''
				holder = []
				for x in line.split():
					holder.append(int(x))
					holder.append(int(x))
					towrite = towrite + x + ' '
				new_times.append([x for x in holder])
				outfile.write(towrite+towrite+'\n')
	outfile.close()
	print 'Got',new_times
	sys.stdout.flush()

	#Check the halt file
	halt = CheckHaltFile(haltfilename)

	#Check if any progress was made
	if halt == 0:
		if old_times == new_times:
			print 'No progress made. Sleeping...'
			sys.stdout.flush()
			time.sleep(10*60)
			halt = CheckHaltFile(haltfilename)
		else:
			print 'Going for the next round'
			sys.stdout.flush()
			old_times = new_times

print 'Halt signal received'

