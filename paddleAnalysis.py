import numpy as np


# This is an analysis script to identify the fraction of events
# with paddle conincidences in the data file. This was motivated by the fact that
# the SiPMs installed by the Krakow group were showing really low rates
# in xh


# import file
# calculate quantities such as NoEvents, total file time, and the average rate
filename = 'cosmicsWithAttenuator2018.dat'
content = open(filename).readlines()
NoEvents = int((len(content)-25)/26)
fileTimeMin = (int(content[-26].split()[1])-int(content[25].split()[1]))/60
aveRate = NoEvents/(fileTimeMin*60)





# this function stores all the paddle information present in an event from the file
# the output is a list with two nested lists. The first list
# has the numerical identifier of the outer paddle(s) that were hit in the event
# the outer paddle identifiers are 1, 2, , ... , 8 for a total of 8 big paddles
# the other list contains the numerical identifier of the inner paddles that were
# hit in the event. the inner paddle identifiers are  1, 2, , ... , 16 for a
# total of 16 small paddles

def getAllPaddleHits(eventNo):
	outerPaddles = np.array(content[36+26*eventNo].split()).astype(int)[1:-1]<1200
	innerPaddles1 = np.array(content[38+26*eventNo].split()).astype(int)[1:-1]<1200
	innerPaddles2 = np.array(content[40+26*eventNo].split()).astype(int)[1:-1]<1200
	innerPaddlesNo = np.concatenate((np.where(innerPaddles1)[0]+1,
		np.where(innerPaddles2)[0]+9))
	return [np.where(outerPaddles)[0]+1,innerPaddlesNo]



# evaluate the function getAllPaddleHits() for all the events in the file
paddles=list(map(getAllPaddleHits,np.arange(0,NoEvents)))



# this function checks which events have two big paddle hits
def twoBigPaddleHitsCheck(eventNo):
	return len(paddles[eventNo][0]) == 2


# get the total number of counts with two big paddle hits
twoBigPaddleCounts =  np.sum(list(map(twoBigPaddleHitsCheck,np.arange(0,NoEvents))))




# this function checks which events have one big paddle hit
def oneBigPaddleHitCheck(eventNo):
	return len(paddles[eventNo][0]) == 1

# get the total number of counts with one big paddle hit
oneBigPaddleCounts =  np.sum(list(map(oneBigPaddleHitCheck,np.arange(0,NoEvents))))



# this function checks which events have a single paddle count
# the single paddle count is defined to be a coincidence between a big paddle
# and one of the two small paddles right in front of the big paddle
def singlePaddleCheck(eventNo):
	if len(paddles[eventNo][0]) == 1:
		outerPaddle1 = paddles[eventNo][0][0]
		return (outerPaddle1*2 in paddles[eventNo][1] or outerPaddle1*2-1 in paddles[eventNo][1])
	else:
		return False



# this function checks which events have at least a single paddle count
# there may be more paddle hits in the events in addition to the
# single paddle count
def singlePaddleCheckOrMore(eventNo):
	if len(paddles[eventNo][0]) >= 1:
		outerPaddle1 = paddles[eventNo][0][0]
		return (outerPaddle1*2 in paddles[eventNo][1] or outerPaddle1*2-1 in paddles[eventNo][1])
	else:
		return False


# get the total number of counts with at least a single paddle
singlePaddleCountsOrMore =  np.sum(list(map(singlePaddleCheckOrMore,np.arange(0,NoEvents))))


# get the total number of counts with a single paddle
singlePaddleCounts =  np.sum(list(map(singlePaddleCheck,np.arange(0,NoEvents))))



# this function checks which events have a double paddle count
# the double paddle count is defined to be a coincidence of two single
# paddle counts
def doublePaddleCheck(eventNo):
	if len(paddles[eventNo][0]) == 2:
		outerPaddle1 = paddles[eventNo][0][0]
		outerPaddle2 = paddles[eventNo][0][1]
		return (outerPaddle1*2 in paddles[eventNo][1] or outerPaddle1*2-1 in paddles[eventNo][1]) and (outerPaddle2*2 in paddles[eventNo][1] or outerPaddle2*2-1 in paddles[eventNo][1])
	else:
		return False


# get the total number of counts with double paddle 
doublePaddleCounts =  np.sum(list(map(doublePaddleCheck,np.arange(0,NoEvents))))


# this function returns True if the event satisfied the new fiber trigger
# condition
def getNewFiberTriggerCondition(eventNo):
	layer1Counts = int(content[41+26*eventNo].split()[1])
	layer2Counts = int(content[43+26*eventNo].split()[1])
	layer3Counts = int(content[45+26*eventNo].split()[1])
	layer4Counts = int(content[48+26*eventNo].split()[1])
	return (layer1Counts or layer2Counts) >= 1 and (layer3Counts or layer4Counts) >= 1


# get all the events in the file that satisfy the new fiber trigger conditon
fiberTriggerCounts=np.sum(list(map(getNewFiberTriggerCondition,np.arange(0,NoEvents))))



# get precise time of event in seconds
def getTime(eventNo):
	return (int(content[25+26*eventNo].split()[1])+
		int(content[25+26*eventNo].split()[2])*10**-6)



# bin the events for histogram
def binEvents(noOfBins):
	arr = np.zeros((noOfBins,2))
	arr[:,0] = np.linspace(0,getTime(NoEvents-1)-getTime(0),noOfBins)
	for eventNo in range(NoEvents):
		arr[int((getTime(eventNo)-getTime(0))/arr[1,0]),1] += 1
	return arr[:-1]


# write fractions and results to file
with open('output.out','a') as f:
	print('file name:', filename,
		'\n','number of events:',NoEvents,
		'\n', 'fraction of events with doublePaddle:',format(doublePaddleCounts/NoEvents,'.5f'),
		'\n', 'fraction of events with twoBigPaddleHits:',format(twoBigPaddleCounts/NoEvents,'.5f'),
		'\n', 'fraction of events with oneBigPaddleHits:',format(oneBigPaddleCounts/NoEvents,'.5f'),
		'\n', 'fraction of events with singlePaddle:',format(singlePaddleCounts/NoEvents,'.5f'),
		'\n', 'fraction of events with at least a singlePaddle:',format(singlePaddleCountsOrMore/NoEvents,'.5f'),
		'\n', 'fraction of events satisfying new fiber trigger:',format(fiberTriggerCounts/NoEvents,'.5f'),
		'\n', 'fileTime in minutes:', format(fileTimeMin, '.1f'),
		'\n', 'Average Rate (Hz) in entire file:',format(aveRate, '.1f'),
		'\n',
		'\n',file=f)




	