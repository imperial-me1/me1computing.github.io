import numpy as np, matplotlib.pyplot as plt, json, os
from stl import mesh
from mpl_toolkits import mplot3d

fileName = 'armReduced.stl'
useJSON = True

def loadMesh(fileName):
    if fileName[-4:].lower() != '.stl':
        fileName += '.stl'
    return mesh.Mesh.from_file(fileName)

def plotPolyMesh(meshObject):
    figure = plt.figure()
    axes = mplot3d.Axes3D(figure)

    # Load the STL files and add the vectors to the plot
    axes.add_collection3d(mplot3d.art3d.Poly3DCollection(armMesh.vectors))

    # Auto scale to the mesh size
    scale = meshObject.points.flatten()
    axes.auto_scale_xyz(scale, scale, scale)

    # Show the plot to the screen
    #plt.show()
    return

def mapVectors(meshObject):

    pointArray, relDict = [], {}

    print('Isolating vectors and building relationship groups:')
    vectorList = meshObject.vectors.tolist()
    printPoints = [int(i) for i in np.linspace(0, len(vectorList)-1,11)]
    for faceNumber in range(0, len(vectorList),1):
        face = vectorList[faceNumber]
        localIndexes = [None,None,None]
        for i in range(0,3,1):
            vector = face[i]
            if not vector in pointArray:
                pointArray.append(vector)
                localIndexes[i] = len(pointArray) - 1
            else:
                localIndexes[i] = pointArray.index(vector)
            continue

        for i in range(0,3,1):
            otherVectorsInFace = localIndexes[:i] + localIndexes[i + 1:]
            try:
                relDict[str(localIndexes[i])].update(otherVectorsInFace)
            except:
                relDict[str(localIndexes[i])] = set(otherVectorsInFace)
            continue
        if faceNumber in printPoints:
            progressIndex = printPoints.index(faceNumber)
            print("["+("="*progressIndex)+("-"*(10-progressIndex))+"]" + " " + str(progressIndex) + "0%")
        continue

    print('Reformating relationship groups')
    for i in range(0,len(pointArray), 1):
        relDict[str(i)] = list(relDict[str(i)])
        continue

    pointArray = np.array(pointArray)
    return pointArray, relDict

def zeroVectors(pointArray):
    print('Zeroing vectors')
    for i in range(0,3,1):
        pointArray[:,i] -= pointArray[:,i].min()
        continue
    return pointArray

def manageMesh(fileName):
    global useJSON
    if os.path.exists(f"mappedVectors_{fileName}.json") and useJSON:
        with open(f"./mappedVectors_{fileName}.json", "r") as file:
            jsonObject = json.load(file)
        pointArray, relDict = np.array(jsonObject["pointArray"]), jsonObject["relDict"]
    else:
        print('Generating vector map:')
        meshObject = loadMesh(fileName)
        pointArray, relDict = mapVectors(meshObject)
        pointArray = zeroVectors(pointArray)
        if useJSON:
            print('Writing to .json')
            with open(f"./mappedVectors_{fileName}.json", "w") as file:
                json.dump({"pointArray": pointArray.tolist(), "relDict": relDict}, file)
    return pointArray, relDict

def getAllRadii(pointArray, relDict):
    numPoints, radiiList, invalidPoints = len(pointArray), [], []
    printPoints = [int(i) for i in np.linspace(0, numPoints-1,11)]
    print('Getting Radii:')
    for i in range(0,numPoints,1):
        relList = relDict[str(i)]
        radiiList.append([])
        xi = pointArray[i]
        for j in range(0,len(relList)-1,1):
            xj = pointArray[relList[j]]
            for k in relList[j+1:]:
                xk = pointArray[k]

                ximxj , ximxk = xi - xj , xi - xk
                angle = np.arccos(np.dot(ximxj / np.linalg.norm(ximxj) , ximxk / np.linalg.norm(ximxk) ) )
                if angle > 2.513:
                    c = 2 * np.sin( angle ) / np.linalg.norm(xj - xk)
                    if c == 0:
                        print("c = 0 for one combination in point", i)
                        radiiList[-1].append(999999)
                    else:
                        radiiList[-1].append(1/c)
                else:
                    radiiList[-1].append(999999)

                continue
            continue

        if len(radiiList[-1]) > 0:
            workingMin = min(radiiList[-1])

            radiiList[-1] = workingMin
            if workingMin == 999999:
                invalidPoints.append(i)
            else:
                pass
        else:
            print("FLAG on point", i)

        if i in printPoints:
            progressIndex = printPoints.index(i)
            print("["+("="*progressIndex)+("-"*(10-progressIndex))+"]" + " " + str(progressIndex) + "0%")
        continue
    print('Filtering radii')
    remainingIndices = np.array([str(i) for i in range(0, numPoints,1)])
    radiiList = np.array(radiiList)
    pointArray, radiiList, remainingIndices = np.delete(pointArray, np.array(invalidPoints), 0), np.delete(radiiList, np.array(invalidPoints), 0), np.delete(remainingIndices, np.array(invalidPoints))
    m = 2
    d = np.abs(radiiList - np.median(radiiList))
    mdev = np.median(d)
    s = d/mdev if mdev else 0.
    return pointArray[s<m], radiiList[s<m], remainingIndices[s<m]

def smoothRadii(radiiList, remainingIndices, relDict):
    smoothedRadii = []
    remainingIndices = remainingIndices.tolist()
    numRadii = len(radiiList)
    printPoints = [int(i) for i in np.linspace(0, numRadii-1, 11)]
    print('Smoothing radii:')
    for i in range(0,numRadii,1):
        strPointNum = remainingIndices[i]
        localPoints = relDict[strPointNum]
        numValidLocalPoints, sumValidLocalPoints = 0, 0
        for j in localPoints:
            if str(j) in remainingIndices:
                numValidLocalPoints += 1
                sumValidLocalPoints += radiiList[int(remainingIndices.index(str(j)))]
            else:
                pass
            continue
        smoothedRadii.append(((radiiList[i]*numValidLocalPoints) + sumValidLocalPoints)/(numValidLocalPoints+1))
        if i in printPoints:
            progressIndex = printPoints.index(i)
            print("["+("="*progressIndex)+("-"*(10-progressIndex))+"]" + " " + str(progressIndex) + "0%")
        continue
    return smoothedRadii

def writeRadiiToFile(arrayToWrite):
    with open(f"./radiiDump.txt", "w") as file:
        file.write("\n".join(str(item) for item in arrayToWrite.tolist()))
    return

def plotPlotCloud(pointArray, smoothedRadii):
    ax = plt.axes(projection='3d')
    ax.scatter(pointArray[:,0], pointArray[:,1], pointArray[:,2], c=smoothedRadii, cmap='plasma', linewidth=0.5)
    scale = [pointArray.min(), pointArray.max()]
    ax.auto_scale_xyz(scale, scale, scale)
    return

def main():
    global fileName

    pointArray, relDict = manageMesh(fileName)
    pointArray, radiiList, remainingIndices = getAllRadii(pointArray, relDict)

    writeRadiiToFile(radiiList)

    smoothedRadii = smoothRadii(radiiList, remainingIndices, relDict)

    plotPlotCloud(pointArray, smoothedRadii)
    plt.show()
    return

main()
