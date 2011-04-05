def stripNumbers(inputstr, method):
	length = len(inputstr)
	strcounter = 0
	tempstr = ""
	returnlist = []
	while strcounter < length:
		if (inputstr[strcounter] < '0' or inputstr[strcounter] > '9') and inputstr[strcounter] != '.':
			tempstr = ""
		else:
			tempstr = tempstr + inputstr[strcounter]
			if strcounter + 1 < length:
				if (inputstr[strcounter + 1] < '0' \
					or inputstr[strcounter + 1] > '9') and inputstr[strcounter + 1] != '.':
						if method == 0:
							returnlist.append(float(tempstr))
						if method == 1:
							returnlist.append(int(tempstr))
						if method == 2:
							returnlist.append(tempstr)
			if strcounter + 1 == length:
					if method == 0:
						returnlist.append(float(tempstr))
					if method == 1:
						returnlist.append(int(tempstr))
					if method == 2:
						returnlist.append(tempstr)
		strcounter = strcounter + 1;	
	return returnlist

