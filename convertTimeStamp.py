
def convertTimeStamp(timeToConvert): 
       
    timeCur = str(timeToConvert)
    sz =len(timeCur)
        
    timeStr = ''
    for i in range(sz,3):
        timeStr+= '0'
    
    timeStr+=timeCur
   
    return(timeStr)



