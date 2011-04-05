def compare(num1, num2, numsigdig):
	tempnum1 = abs(num1)
	tempnum2 = abs(num2)
	c = 0
	index = 0
	reachdigit = 0
	returnval = False
	if num1 < 0 and num2 > 0:
		return False
	if num1 > 0 and num2 < 0:
		return False
	if num1 == num2:
		return True
	if tempnum1 < 1 and tempnum2 > 1:
		return False
	if tempnum1 > 1 and tempnum2 < 1:
		return False
	tempnum1 = str(tempnum1)
	tempnum2 = str(tempnum2)
	len1 = len(tempnum1)
	len2 = len(tempnum2)
	while c < numsigdig:
		if index < len1 and index < len2:

			if tempnum1[index] == '.' and tempnum2[index] == '.':
				index = index + 1
				continue
			if tempnum1[index] == tempnum2[index] and reachdigit == 0 and (tempnum1[index] == '0' or tempnum1[index] == '.'):
				index = index + 1
				continue
			if tempnum1[index] == tempnum2[index] and reachdigit == 0 and tempnum1[index] != '0' and \
													tempnum1[index] != '.':
				reachdigit = 1
				index = index + 1
				c = c + 1
				continue
			if tempnum1[index] == tempnum2[index] and reachdigit == 1:
				index = index + 1
				c = c + 1
				continue
			if tempnum1[index] != tempnum2[index]:
				return False
		else:
			print "Warning : the precision of input data is not enough"
			return False
	return True
		
