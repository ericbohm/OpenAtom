def stripNumbers(inputstr, method):
	length = len(inputstr)
	strcounter = 0
	tempstr = ""
	returnlist = []
#start to check every char inside the input string
	while strcounter < length:
# if we still haven't reach the a digit then simply empty the tempstr
		if (inputstr[strcounter] < '0' or inputstr[strcounter] > '9') and inputstr[strcounter] != '.':
			tempstr = ""
		else:
# else, store the digit into the tempstr
			tempstr = tempstr + inputstr[strcounter]
			if strcounter + 1 < length:
# if we haven't reach the end of the string and the next char is not a digit, which means we have already strip a number from the string into tempstr, then we need to store it in to the list.
				if (inputstr[strcounter + 1] < '0' or inputstr[strcounter + 1] > '9') and inputstr[strcounter + 1] != '.':
						if method == 0:
							returnlist.append(float(tempstr))
						if method == 1:
							returnlist.append(int(tempstr))
						if method == 2:
							returnlist.append(tempstr)
# if we reach the end of the string, which means we have already strip a number from the string into tempstr then we need to store it in to the list.
			if strcounter + 1 == length:
					if method == 0:
						returnlist.append(float(tempstr))
					if method == 1:
						returnlist.append(int(tempstr))
					if method == 2:
						returnlist.append(tempstr)
# if we haven't reach the end of the string and the next char is also a digit then simply do next iteration
		strcounter = strcounter + 1;	
	return returnlist

