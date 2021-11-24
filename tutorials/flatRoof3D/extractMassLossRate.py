import numpy as np
import matplotlib.pyplot as plt

def extractData(nProc, time, dataName):
    dataList = []
    for i in range(nProc):
        folderName = 'processor' + str(i)
        with open(folderName + '/' + str(time) + '/' + dataName, 'r') as f:
            flagIndex = False
            dataCount = 0
            for line in f:
                if not flagIndex:
                    try:
                       num = int(line[:-1])
                       #print("%d faces are found in processor %d."%(num, i))
                       flagIndex = True
                    except:
                        pass
                    continue
                try:
                    dataList.append(float(line[:-1]))
                    dataCount += 1
                except:
                    pass

            #print("%d datas are found in processor %d.\n"%(dataCount, i))
    return np.array(dataList)
    
def main():
    timeSeries = [1000*i for i in range(1000)]
    massLoss = []
    count= 0
    with open('mass.csv', 'w') as f:
        f.write("Time[s], Mass[kg]\n")
        firstStep = None
        for i in timeSeries[1:]:
            try:
                deltaH = extractData(16, i, 'CellCentreZ')
                areaZ = extractData(16, i, 'zArea')
            except FileNotFoundError:
                continue
            if firstStep is None:
                firstStep = i
            count += 1
            massLoss.append(np.sum(areaZ * (deltaH - 8)))
            f.write("%d, %f\n"%(i, massLoss[-1]*150))
    #massLoss = np.array(massLoss)
    #plt.plot(timeSeries[1:], massLoss*150)
    #plt.plot(timeSeries[2:count+1], massLoss[:-1] - massLoss[1:])
    #plt.show()
if __name__ == '__main__':
    main()
