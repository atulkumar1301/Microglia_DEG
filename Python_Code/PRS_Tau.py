#This code is a pipeline for Cell-Specific PRS calculation for continuous Variable
import subprocess
import os

Cell_Type = ["Original_PRS"]
plink2 = "/Volumes/ATUL_6TB/Tools/./plink2"
plink = "/Volumes/ATUL_6TB/Tools/plink_mac_20230116/./plink"
base_file = "/Volumes/ATUL_6TB/Work/Projects/Microglia_DEG/BF_2/GWAS_Microglia_BF-2"
GWAS_Summary = "/Volumes/ATUL_6TB/Data/GWAS_Summary_Statistics/AD_GWAS_Bellenguez_2022.txt"
pwd = "/Volumes/ATUL_6TB/Work/Projects/Microglia_DEG/BF_2/"
f_n = "5_Data_Full_Imputed_Analysis_Taupet.txt"

##Extracting the Sample ID
print ("Id Extraction Started")
Out_File = open ("/Volumes/ATUL_6TB/Work/Projects/Microglia_DEG/BF_2/Keep_Patient.txt", 'w', 1)
with open ("/Volumes/ATUL_6TB/Work/Projects/Microglia_DEG/BF_2/" + f_n, 'r') as datafile:
    IDs = datafile.readline ()
    IDs_Split = IDs.split ("\t")
    Out_File.write (str (IDs_Split [0]) + "\t" + str (IDs_Split [1]) + "\n")
    for IDs in datafile:
        IDs_Split = IDs.split ("\t")
        Out_File.write (str (IDs_Split [0]) + "\t" + str (IDs_Split [1]) + "\n")
print ("Id Extraction finished")


print ("Extracting GWAS data for Sample")
subprocess.run ([plink2, "--pfile", base_file, "--keep", "Keep_Patient.txt", "--maf", "0.05", "--allow-extra-chr", "--make-bed", "-out", "9_QC_GWAS_data"], cwd = pwd)
print ("Extraction of GWAS data for Sample Complete")

print ("Calculating the PCA of genetic Data")
subprocess.run ([plink2, "--bfile", "9_QC_GWAS_data", "--allow-extra-chr", "--pca", "--out", "PCA_FILE"], cwd = pwd)
print ("Calculation of the PCA of genetic Data Complete")

print ("Started the Clumping")
subprocess.run ([plink, "--bfile", "9_QC_GWAS_data", "--allow-extra-chr", "--clump", GWAS_Summary, "--clump-kb", "1000", "--clump-p1", "1", "--clump-p2", "1", "--clump-r2", "0.1", "--out", "Clumped_File"], cwd = pwd)
print ("Finished Clumping")

print ("Calling R for making Pre PRS File")
subprocess.call ("/Volumes/ATUL_6TB/Tools/R_Codes/PRS.R", cwd = pwd)
print ("Generation of Pre PRS File Complete")

for cell in Cell_Type:
    if cell == "Original_PRS":
        print ("Started Generating Score file")
        path = pwd + cell + "/"
        C_Main = os.path.exists (path)
        if not C_Main:
            os.mkdir (path)
        APOE = path + "APOE"
        C_1 = os.path.exists (APOE)
        if not C_1:
            os.mkdir (APOE)
        Non_APOE = path + "Non_APOE"
        C_2 = os.path.exists (Non_APOE)
        if not C_2:
            os.mkdir (Non_APOE)
        f_m = open (path + "SNP_Extract.txt", 'w', 1)
        with open (pwd + "Pre_PRS.txt", 'r') as myFile:
            line = myFile.readline ()
            for line in myFile:
                line_list = line.split("\t")
                f_m.write (str (line_list [2]) + "\n")
        SNP_Extract = path + "SNP_Extract.txt"
        out_file = path + "10_QC_GWAS_data"
        subprocess.run ([plink2, "--bfile", "9_QC_GWAS_data", "--extract", SNP_Extract, "--allow-extra-chr", "--make-pgen", "--out", out_file], cwd = pwd)
        l = [0.05, 0.005, 0.0005, 0.00005, 0.000005, 0.0000005, 0.00000005]
        f_m_1 = open(path + "APOE_Region_SNPs.txt", 'w', 1)
        with open (pwd + "Pre_PRS.txt", 'r') as myFile_1:
            line_1 = myFile_1.readline ()
            f_m_1.write (line_1)
            for line_1 in myFile_1:
                line_list_1 = line_1.split ("\t")
                if int (line_list_1 [0]) == 19:
                    if ((int (line_list_1 [1]) > 43905790) and (int (line_list_1 [1]) < 45909394)):
                        f_m_1.write (line_1)

### Including APOE Region Variant
        for i in l:
            f_m_2 = open(APOE + "/p_value_" + str (i) + ".txt", 'w', 1)
            with open (pwd + "Pre_PRS.txt") as myFile_2:
                line_2 = myFile_2.readline ()
                for line_2 in myFile_2:
                    line_list_2 = line_2.split ("\t")
                    if float (line_list_2 [5]) <= i:
                        f_m_2.write (line_list_2 [2] + "\t" + str(line_list_2 [3]).upper()+ "\t" + line_list_2 [4] + "\n")
        for i in l:
            f_i = APOE + "/p_value_" + str(i) + ".txt"
            f_o = APOE + "/p_value_" + str(i)
            subprocess.run ([plink2, "--pfile", "10_QC_GWAS_data", "--score", f_i, "--allow-extra-chr", "--out", f_o], cwd = path)
        R_code = subprocess.call (["/Volumes/ATUL_6TB/Tools/R_Codes/PRS_Linear.R", str(f_n)], cwd = APOE)

### Excluding APOE Region Variant
        for i in l:
            f_m_3 = open(Non_APOE + "/p_value_" + str (i) + ".txt", 'w', 1)
            with open (pwd + "Pre_PRS.txt") as myFile_3:
                line_3 = myFile_3.readline ()
                for line_3 in myFile_3:
                    line_list_3 = line_3.split ("\t")
                    if int (line_list_3 [0]) == 19:
                        if ((int (line_list_3 [1]) > 43905790) and (int (line_list_3 [1]) < 45909394)):
                            print ("PASS")
                        elif float (line_list_3 [5]) <= i:
                            f_m_3.write (line_list_3 [2] + "\t" + line_list_3 [3] + "\t" + line_list_3 [4] + "\n")
                    elif float (line_list_3 [5]) <= i:
                        f_m_3.write (line_list_3 [2] + "\t" + str(line_list_3 [3]).upper()+ "\t" + line_list_3 [4] + "\n")
        for i in l:
            f_i = Non_APOE+ "/p_value_" + str(i) + ".txt"
            f_o = Non_APOE + "/p_value_" + str(i)
            subprocess.run ([plink2, "--pfile", "10_QC_GWAS_data", "--score", f_i, "--allow-extra-chr", "--out", f_o], cwd = path)
        R_code = subprocess.call (["/Volumes/ATUL_6TB/Tools/R_Codes/PRS_Linear_Non_APOE.R", str(f_n)], cwd = Non_APOE)
        print ("Finished Generating Score file")

