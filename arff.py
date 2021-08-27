from __future__ import print_function
from copy import deepcopy
import sys
import random
import math
import time
import datetime
import os
import subprocess
import collections

Amino_Acid_1 = ('G','A','V','L','I','F','W',
	'Y','D','H','N','E','K','Q',
	'M','R','S','T','C','P')

Amino_Acid_2 = deepcopy(Amino_Acid_1)

Another_Amino_Acid = 'X'

Hydrophobicity = [['R','K','E','D','Q','N'],['G','A','S','T','P','H','Y'],['C','L','V','I','M','F','W']]

Normalized_van_der_waals_volume = [['G','A','S','T','P','D'], ['N','V','E','Q','I','L'], ['M','H','K','F','R','Y','W']]

Polarity = [['L','I','F','W','C','M','V','Y'], ['P','A','T','G','S'], ['H','Q','R','K','N','E','D']]

Polarizability = [['G','A','S','D','T'], ['C','P','N','V','E','Q','I','L'], ['K','M','H','F','R','Y','W']]

Charge = [['K','R'], ['A','N','C','Q','G','H','I','L','M','F','P','S','T','W','Y','V'], ['D','E']]

Secondary_structure = [['E','A','L','M','Q','K','R','H'], ['V','I','Y','C','W','F','T'], ['G','N','P','S','D']]

Solvent_accessibility = [['A','L','F','C','G','I','V','W'], ['R','K','Q','E','N','D'], ['M','P','S','T','H','Y']]

Properties_List = [ Hydrophobicity, Normalized_van_der_waals_volume, Polarity, Polarizability, Charge, Secondary_structure, Solvent_accessibility ]

Properties_String = [ 'Hydrophobicity', 'Normalized_van_der_waals_volume', 'Polarity', 'Polarizability', 'Charge', 'Secondary_structure', 'Solvent_accessibility' ]

Front_to_End_String = ['_Front_To_End_10_','_Front_To_End_20_','_Front_To_End_30_','_Front_To_End_40_','_Front_To_End_50_','_Front_To_End_60_','_Front_To_End_70_','_Front_To_End_80_','_Front_To_End_90_','_Front_To_End_100_']

Class_String = [ 'Class_1', 'Class_2', 'Class_3' ]



def Phy_sci_Statistics(seq,Attribute_List,Attribute_String):

	Statistics = dict()
	Sequence_Quantity = 0
	for Amino_Acid_Name,Amino_Acid_Sequence in seq.iteritems():
		Dict_Key_Pre_String = '_Front_To_End_'
		Sequence_To_Handle = Amino_Acid_Sequence
		Float_Front_to_end_Calculate = 0
		Key_Front_end_Count = 0			#Count the 0.1 ~ 1.0.
		Float_Parameter = 0.1			#
		Key_Parameter = 10
		# Final_Turn = True
		# Reverse_Turn = True
		Sequence_Quantity += 1
		Initial_String_Length = 0		#Calculate the start Amino Acid.
		Partial_String_Length = 0		#Calculate the end of Amino Acid.
		while(True):
			Class_1 = 0
			Class_2 = 0
			Class_3 = 0
			Float_Front_to_end_Calculate += Float_Parameter
			Key_Front_end_Count += Key_Parameter
			# print ( str('%.1f'%Front_to_end) + ' ' + Dict_Key)
			# print ( round(10 / Front_to_end) )

			Length = len(Sequence_To_Handle)
			Initial_String_Length = Partial_String_Length
			Partial_String_Length = int(round(Length * Float_Front_to_end_Calculate))
			# print (Float_Front_to_end_Calculate)

			# Reverse======================================================
			# if not Final_Turn and Reverse_Turn:
			# 	Sequence_To_Handle = Sequence_To_Handle[::-1]
			# 	Dict_Key_Pre_String = 'Hydrophobicity_End_To_Front_'
			# 	Reverse_Turn = False

			Dict_Key_C1 = Attribute_String + Dict_Key_Pre_String + str(Key_Front_end_Count) + '_Class_1'
			Dict_Key_C2 = Attribute_String + Dict_Key_Pre_String + str(Key_Front_end_Count) + '_Class_2'
			Dict_Key_C3 = Attribute_String + Dict_Key_Pre_String + str(Key_Front_end_Count) + '_Class_3'

			Partial_String = Sequence_To_Handle[Initial_String_Length:Partial_String_Length]
			# print("Before: " + Amino_Acid_Sequence)
			# print ("After: " + Sequence_To_Handle)
			# print("Partial_" + str(Initial_String_Length) + "~" + str(Partial_String_Length) + ": " + Partial_String)

			for Amino_Acid in Partial_String:
				if Amino_Acid in Attribute_List[0]:
					Class_1 += 1
				elif Amino_Acid in Attribute_List[1]:
					Class_2 +=1
				elif Amino_Acid in Attribute_List[2]:
					Class_3 +=1
			if Statistics.has_key(Dict_Key_C1) and Statistics.has_key(Dict_Key_C2) and Statistics.has_key(Dict_Key_C3):
				Statistics[Dict_Key_C1] += float('%.3f'%(Class_1 / float(Length)))
				Statistics[Dict_Key_C2] += float('%.3f'%(Class_2 / float(Length)))
				Statistics[Dict_Key_C3] += float('%.3f'%(Class_3 / float(Length)))
			else:	
				Statistics[Dict_Key_C1] = float('%.3f'%(Class_1 / float(Length)))
				Statistics[Dict_Key_C2] = float('%.3f'%(Class_2 / float(Length)))
				Statistics[Dict_Key_C3] = float('%.3f'%(Class_3 / float(Length)))
			# Debug========================================================================
			# print(str(Hydrophobicity[0]) + ": " + str(Class_1) + ", " + '%.2f'%(Class_1 / float(Length)) + ", " + str(Length))
			# print(str(Hydrophobicity[1]) + ": " + str(Class_2) + ", " + '%.2f'%(Class_2 / float(Length)) + ", " + str(Length))
			# print(str(Hydrophobicity[2]) + ": " + str(Class_3) + ", " + '%.2f'%(Class_3 / float(Length)) + ", " + str(Length))



			# if str('%.1f'%Float_Front_to_end_Calculate) == '1.0':
			# 	Final_Turn = False
			# 	Float_Parameter = -0.1
			# 	Key_Parameter = -10
			# if not Final_Turn and str('%.1f'%Float_Front_to_end_Calculate) == '0.1':
			# 	break

			if str('%.1f'%Float_Front_to_end_Calculate) == '1.0':
				break

	for Dict_Key,Data in Statistics.iteritems():
		Temp = (Data / Sequence_Quantity)
		# Statistics[Dict_Key] = '%.2f'%math.log(Temp,2)
		Statistics[Dict_Key] = '%.3f'%Temp

	return Statistics

def AAC_AAPC_statistics(seq):
	
	AAC_AAPC_Statistics = dict()
	Total_number = 0
	for Amino_Acid_Name,Amino_Acid_Sequence in seq.iteritems():
		Total_number += 1
		for AA_1 in Amino_Acid_1:
			for AA_2 in Amino_Acid_2:
				AAP = AA_1 + AA_2

			# =============================AAPC Statisitcs.=============================#
				if AAP in Amino_Acid_Sequence:											# AAP belong to Amino Acid Sequence.
					AAPnum = 0															#
					for i in range(len(Amino_Acid_Sequence) - 1):						# Calculate the AAP count.
						temp = Amino_Acid_Sequence[i]+Amino_Acid_Sequence[i+1]			#
						if(temp == AAP):												# Which in the Animo Acid Sequence.
							AAPnum += 1													#
					result = float("%.4f"% (AAPnum/float(len(Amino_Acid_Sequence) - 1)))
					
					if AAP in AAC_AAPC_Statistics:										# AAPC belong to Amino Acid Sequence.
						AAC_AAPC_Statistics[AAP] = AAC_AAPC_Statistics[AAP] + result	# AAPC has exist in AAC_AAPC_Statistics.

					else:																# AApc belong to Amino Acid Sequence.
						AAC_AAPC_Statistics[AAP] = result								# AAPC doesn't exist in AAC_AAPC_Statistics yet.
				
			#	else:																	# AAP doesn't belong to Amino Acid Sequence.
			#		AAC_AAPC_Statistics[AAP] = 0.0										#
			# =========================== End AAPC Statisitcs.==========================#


			# ==============================AAC Statisitcs.=============================#

			if AA_1 in Amino_Acid_Sequence:												# AA_1 belong to Amino Acid Sequence.
				result = float("%.4f"% (int(Amino_Acid_Sequence.count(AA_1))/float(len(Amino_Acid_Sequence))))
				
				if AA_1 in AAC_AAPC_Statistics:											# AA_1 belong to Amino Acid Sequence.
					AAC_AAPC_Statistics[AA_1] = AAC_AAPC_Statistics[AA_1] + result		# And AA_1 has exist in AAC AAPC Statistics.

				else:																	# AA_1 belong to Amino Acid Sequence.
					AAC_AAPC_Statistics[AA_1] = result									# But AA_1 doesn't exist in AAC AAPC Statistics.

			#else:																		# AA_1 doesn't belong to Amino Acid Sequence.
			#	AAC_AAPC_Statistics[AA_1] = 0.0											# And then I will get AAC_AAPC_Statistics
																						# Which include all of Amino Acid and Amino Acid pair.
			# ============================ End AAC Statisitcs.==========================#
		Another_Amino_Acid_Count = 0
		for AA_X in Amino_Acid_Sequence:
			if AA_X not in Amino_Acid_1:
				Another_Amino_Acid_Count +=1
		Another_Amino_Acid_Result = (Another_Amino_Acid_Count / float(len(Amino_Acid_Sequence)))
		if Another_Amino_Acid in AAC_AAPC_Statistics:											# AA_1 belong to Amino Acid Sequence.
			AAC_AAPC_Statistics[Another_Amino_Acid] = AAC_AAPC_Statistics[Another_Amino_Acid] + Another_Amino_Acid_Result		# And AA_1 has exist in AAC AAPC Statistics.
		else:																	# AA_1 belong to Amino Acid Sequence.
			AAC_AAPC_Statistics[Another_Amino_Acid] = Another_Amino_Acid_Result

	for Amino_Acid_Name in AAC_AAPC_Statistics:
		Temp = AAC_AAPC_Statistics[Amino_Acid_Name] / Total_number
		if Amino_Acid_Name in Amino_Acid_1 :
			AAC_AAPC_Statistics[Amino_Acid_Name] = str("%.3f"%Temp)
		elif Amino_Acid_Name == Another_Amino_Acid:
			AAC_AAPC_Statistics[Amino_Acid_Name] = str("%.5f"%Temp)
		else:
			AAC_AAPC_Statistics[Amino_Acid_Name] = "%.1f"% math.log(Temp,2)
	return AAC_AAPC_Statistics


def Print_ARFF_format(fp):
	fp.write("@relation AMP_Feature_investigation\n")
	fp.write("\n")
	for Properties_Index in range(0,7):
		for mid in Front_to_End_String:
			for rear in Class_String:
				key = Properties_String[Properties_Index] + mid + rear
				fp.write("@attribute " + key + " numeric\n")
	for AA_1 in Amino_Acid_1:											#
		for AA_2 in Amino_Acid_2:										#
			AAP = AA_1+AA_2												# 
			fp.write("@attribute "+AAP+"_Animo_Acid_Pair numeric\n")	# AAPC format.
	
	for AA in Amino_Acid_1:												#
		fp.write("@attribute "+AA+"_Animo_Acid numeric\n")				# AAC format.
	fp.write("@attribute "+Another_Amino_Acid+"_Animo_Acid numeric\n\n")
	fp.write("@attribute class { \"1\" ,\"0\" }\n")
	fp.write("\n")
	fp.write("@data\n\n")


def Phy_sci_ARFF(Amino_Acid_Sequence,Attribute_List,arff_fp,SVM_fp,index):

	Float_Front_to_end_Calculate = 0
	Float_Parameter = 0.1			#

	Initial_String_Length = 0		#Calculate the start Amino Acid.
	Partial_String_Length = 0		#Calculate the end of Amino Acid.
	while(True):
		
		Sequence_To_Handle = Amino_Acid_Sequence
		Class_1 = 0
		Class_2 = 0
		Class_3 = 0
		Float_Front_to_end_Calculate += Float_Parameter
		# print ( str('%.1f'%Front_to_end) + ' ' + Dict_Key)
		# print ( round(10 / Front_to_end) )

		Length = len(Sequence_To_Handle)
		Initial_String_Length = Partial_String_Length
		Partial_String_Length = int(round(Length * Float_Front_to_end_Calculate))
		# print (Float_Front_to_end_Calculate)

		# Reverse======================================================
		# if not Final_Turn and Reverse_Turn:
		# 	Sequence_To_Handle = Sequence_To_Handle[::-1]
		# 	Dict_Key_Pre_String = 'Hydrophobicity_End_To_Front_'
		# 	Reverse_Turn = False
		# Dict_Key_C1 = Attribute_String + Dict_Key_Pre_String + str(Key_Front_end_Count) + '_Class_1'
		# Dict_Key_C2 = Attribute_String + Dict_Key_Pre_String + str(Key_Front_end_Count) + '_Class_2'
		# Dict_Key_C3 = Attribute_String + Dict_Key_Pre_String + str(Key_Front_end_Count) + '_Class_3'
		
		Partial_String = Sequence_To_Handle[Initial_String_Length:Partial_String_Length]
		# print("Before: " + Amino_Acid_Sequence)
		# print ("After: " + Sequence_To_Handle)
		# print("Partial_" + str(Initial_String_Length) + "~" + str(Partial_String_Length) + ": " + Partial_String)

		for Amino_Acid in Partial_String:
			if Amino_Acid in Attribute_List[0]:
				Class_1 += 1
			elif Amino_Acid in Attribute_List[1]:
				Class_2 += 1
			elif Amino_Acid in Attribute_List[2]:
				Class_3 += 1
		result_1 = float('%.3f'%(Class_1 / float(Length)))
		result_2 = float('%.3f'%(Class_2 / float(Length)))
		result_3 = float('%.3f'%(Class_3 / float(Length)))

		arff_fp.write( str(result_1) + ',' )
		index += 1
		SVM_fp.write( str(index)+":"+str(result_1) + ' ')

		arff_fp.write( str(result_2) + ',' )
		index += 1
		SVM_fp.write( str(index)+":"+str(result_2) + ' ')

		arff_fp.write( str(result_3) + ',' )
		index += 1
		SVM_fp.write( str(index)+":"+str(result_3) + ' ')

		# if Statistics.has_key(Dict_Key_C1) and Statistics.has_key(Dict_Key_C2) and Statistics.has_key(Dict_Key_C3):
		# 	Statistics[Dict_Key_C1] += float('%.2f'%(Class_1 / float(Length)))
		# 	Statistics[Dict_Key_C2] += float('%.2f'%(Class_2 / float(Length)))
		# 	Statistics[Dict_Key_C3] += float('%.2f'%(Class_3 / float(Length)))
		# else:	
		# 	Statistics[Dict_Key_C1] = float('%.2f'%(Class_1 / float(Length)))
		# 	Statistics[Dict_Key_C2] = float('%.2f'%(Class_2 / float(Length)))
		# 	Statistics[Dict_Key_C3] = float('%.2f'%(Class_3 / float(Length)))
		# Debug========================================================================
		# print(str(Hydrophobicity[0]) + ": " + str(Class_1) + ", " + '%.2f'%(Class_1 / float(Length)) + ", " + str(Length))
		# print(str(Hydrophobicity[1]) + ": " + str(Class_2) + ", " + '%.2f'%(Class_2 / float(Length)) + ", " + str(Length))
		# print(str(Hydrophobicity[2]) + ": " + str(Class_3) + ", " + '%.2f'%(Class_3 / float(Length)) + ", " + str(Length))

		# if str('%.1f'%Float_Front_to_end_Calculate) == '1.0':
		# 	Final_Turn = False
		# 	Float_Parameter = -0.1
		# 	Key_Parameter = -10
		# if not Final_Turn and str('%.1f'%Float_Front_to_end_Calculate) == '0.1':
		# 	break

		if str('%.1f'%Float_Front_to_end_Calculate) == '1.0':
			break


	return index


def AAC_AAPC_ARFF(arff_fp,SVM_fp,Seq,check,train_or_test):


	for Amino_Acid_Name,Amino_Acid_Sequence in Seq.iteritems():
		
		if(train_or_test):				#
			if(check):					# SVM format
				SVM_fp.write("1" + ' ')	#
			else:						# '1' or '0' represent positive or negative.
				SVM_fp.write("0" + ' ')	# 		
		else:							# '?' represent unknow.
			SVM_fp.write("?" + ' ')		#

		index = 0	
		Catch_Data = dict()
		for Properties_Index in range(0,7):
			index = Phy_sci_ARFF(Amino_Acid_Sequence,Properties_List[Properties_Index],arff_fp,SVM_fp,index)

		# ===========================AAPC Statistics.===========================#
		for AA_1 in Amino_Acid_1:		# AAPC Statistics.
			for AA_2 in Amino_Acid_2:
				AAP = AA_1 + AA_2
				index += 1

				if AAP in Amino_Acid_Sequence:									# AAP belong to the Sequence.
					AAPnum = 0													#
					for i in range(len(Amino_Acid_Sequence) - 1):				# Calculate the AAP count.
						temp = Amino_Acid_Sequence[i]+Amino_Acid_Sequence[i+1]	#
						if(temp == AAP):										# Which in the Animo Acid Sequence.
							AAPnum += 1											#
					result = float("%.4f"% (AAPnum/float(len(Amino_Acid_Sequence) - 1)))
					arff_fp.write( str(result) + ',' ) 							# Arff format
					SVM_fp.write( str(index)+":"+str(result) + ' ')				# SVM format

				else:															# AAP doesn't belong to the Sequence.
					arff_fp.write( "0" + ',')									# Arff format
					SVM_fp.write(str(index)+":"+'0' + ' ')						# SVM format
		# ======================== End AAPC Statistics. ========================#

		# ===========================AAC Statistics.============================#
		for AA in Amino_Acid_1:
			index+=1
			if AA in Amino_Acid_Sequence:										# AA belong to the Sequence.
				result = float("%.4f"% (int(Amino_Acid_Sequence.count(AA))/float(len(Amino_Acid_Sequence))))
				arff_fp.write( str(result) + ',' )								# Arff format
				SVM_fp.write( str(index) + ":" + str(result) + ' ')				# SVM format

			else:																# AA doesn't belong to the Sequence. 
				arff_fp.write( "0" + ',')										# Arff format
				SVM_fp.write(str(index)+":"+'0' + ' ')							# SVM format
		# ========================= End AAC Statistics. ========================#
		Another_Amino_Acid_Count = 0
		for AA_X in Amino_Acid_Sequence:
			if AA_X not in Amino_Acid_1:
				Another_Amino_Acid_Count +=1
		Another_Amino_Acid_Result = float("%.5f"%(Another_Amino_Acid_Count / float(len(Amino_Acid_Sequence))))
	
		index += 1
		arff_fp.write( str(Another_Amino_Acid_Result) + ',' )								# Arff format
		SVM_fp.write( str(index) + ":" + str(Another_Amino_Acid_Result) + ' ')

		SVM_fp.write('\n')

		if(train_or_test):				#
			if(check):					# arff format.
				arff_fp.write("1\n")	#
			else:						# '1' or '0' represent positive or negative.
				arff_fp.write("0\n")	#
		else:							# '?' represent unknow.
			arff_fp.write("?\n")		#

			
###test
localtime = time.localtime(time.time())
#print (localtime.tm_year,localtime.tm_mon,localtime.tm_mday,sep = '_')
#nowtime = str(localtime.tm_year)+"_"+str(localtime.tm_mon)+"_"+str(localtime.tm_mday)+"_"+str(localtime.tm_hour)+str(localtime.tm_min)+str(localtime.tm_sec)
now = datetime.datetime.now()
#nowtime = now.strftime('%Y%m%d')
nowtoday = sys.argv[1]
nowtime = sys.argv[2]

#print(nowtime)
path = "/home/dbamp/public_html/temp/"+nowtoday+"/"+nowtime
#if not os.path.isdir(path):
#	os.makedirs(path)
	
seq = collections.OrderedDict()
tmpName = ''
while(True):
	
	try:
		input = raw_input()
		if(input[0:1] == '>'):	
			tmpName = input[1:].strip()
			#print(tmpName)
		else:
			seq[tmpName] = input.strip()
	except EOFError:
		break;

path1 = "/home/dbamp/public_html/temp/"+nowtoday+"/SVM.txt"
SVM_Train_fp = open(path1,'w')
Arff_Train_fp = open(path+".arff","w")
Print_ARFF_format(Arff_Train_fp)
AAC_AAPC_ARFF(Arff_Train_fp,SVM_Train_fp,seq,True,True)
#AAC_AAPC_ARFF(Arff_Train_fp,SVM_Train_fp,Negative_tr,False,True)
Arff_Train_fp.close()
SVM_Train_fp.close()
time = (nowtoday, nowtime)
#subprocess.call("php /home/AMPrid/public_html/SRPCat/public_html/Amphibia_FS.php nowtime")
os.system("php /home/dbamp/public_html/Amphibia_FS.php %s %s"%(time))
os.system("php /home/dbamp/public_html/Bacteria_FS.php %s %s"%(time))
os.system("php /home/dbamp/public_html/Fish_FS.php %s %s"%(time))
os.system("php /home/dbamp/public_html/Human_FS.php %s %s"%(time))
os.system("php /home/dbamp/public_html/Insects_FS.php %s %s"%(time))
os.system("php /home/dbamp/public_html/Mammals_FS.php %s %s"%(time))
os.system("php /home/dbamp/public_html/Plants_FS.php %s %s"%(time))