set BaseFolder=S:\U_Proteomica\LABS\LAB_ARR\MouseLDL-HighFattyDiet\PTMs\Quant_WF2
set Q2CRelationFile="S:\U_Proteomica\LABS\LAB_ARR\MouseLDL-HighFattyDiet\PTMs\Quant_WF2\protein2category.tsv"
set Data=S:\U_Proteomica\LABS\LAB_ARR\MouseLDL-HighFattyDiet\PTMs\Quant_WF2\ID_q_XV_withProtAccesion.txt
set PDMTable="S:\U_Proteomica\LABS\LAB_ARR\MouseLDL-HighFattyDiet\PTMs\Sticker\V16\ALL_calibrated_DMTable_PeakAssignation_FDR_FDRfiltered_DM0Solved_TrunkSolved_X=6_PeakAssignation_SS_PDMTable_Sticker.txt"

set aljamiaSData=aljamia.exe -x"%Data%" -p"%BaseFolder%\%%i" -o"Scans_uncalibrated.xls" -aS2PDM_inOUTs_uncalibrated -i"[Raw_FirstScan]" -j"[X%%i_Mean]" -k"[V%%i_Mean]" -R1
set aljamiaS2PDMRels=aljamia.exe -x"%Data%" -p"%BaseFolder%\%%i" -o"S2PDM_RelationsFile.xls" -aS2PDM_RelationsFile -i"[SiteSequence]" -j"[Raw_FirstScan]" -R1
set aljamiaPDM2qRels=aljamia.exe -x"%Data%" -p"%BaseFolder%\%%i" -o"PDM2q_RelationsFile.xls" -aPDM2q_RelationsFile -i"[Protein_Accession]" -j"[SiteSequence]" -R1
set aljamiaPDM2qdnaRels=aljamia.exe -x"%PDMTable%" -p"%BaseFolder%\%%i" -o"PDM2qdna_RelationsFile.xls" -aPDM2qdna_RelationsFile -i"[qdna]" -j"[pdm]" -R1
set aljamiaqdna2qnaRels=aljamia.exe -x"%PDMTable%" -p"%BaseFolder%\%%i" -o"qdna2qna_RelationsFile.xls" -aqdna2qna_RelationsFile -i"[q]:[n]:[a]" -j"[qdna]" -R1
set aljamiaqna2qRels=aljamia.exe -x"%PDMTable%" -p"%BaseFolder%\%%i" -o"qna2q_RelationsFile.xls" -aqna2q_RelationsFile -i"[q]" -j"[q]:[n]:[a]" -R1

REM REM REM Basic WF (protein with all PDM)
set klibrate=klibrate.exe -d"%BaseFolder%\%%i\Scans_uncalibrated.xls" -r"%BaseFolder%\%%i\S2PDM_RelationsFile.xls" -p"%BaseFolder%\%%i" -aS2PDM_inOUTs_calibrated -o"scan.xls" -g
set sanxotS2PDM_in_outs=sanxot.exe -aS2PDM_inOuts -p"%BaseFolder%\%%i" -d"scan.xls" -r"%BaseFolder%\%%i\S2PDM_RelationsFile.xls" -g --tags="" --emergencyvariance -V"S2PDM_inOUTs_calibrated_infoFile.txt"
set sanxotsieveSPDM=sanxotsieve.exe -aS2PDMOuts -p"%BaseFolder%\%%i" -d"scan.xls" -r"%BaseFolder%\%%i\S2PDM_RelationsFile.xls" -f0.01 -V"S2PDM_inOuts_infoFile.txt"
set sanxotS2PDM_no_outs=sanxot.exe -aS2PDM_noOuts -p"%BaseFolder%\%%i" -d"scan.xls" -r"%BaseFolder%\%%i\S2PDMOuts_tagged.xls" -o"PDM.xls" -g -V"S2PDM_inOuts_infoFile.txt" -f --tags="!out"

set sanxotPDM2q_in_outs=sanxot.exe -aPDM2q_inOuts -p"%BaseFolder%\%%i" -d"PDM.xls" -r"%BaseFolder%\%%i\PDM2q_RelationsFile.xls" -g --tags=""
set sanxotsievePDMq=sanxotsieve.exe -aPDM2qOuts -p"%BaseFolder%\%%i" -d"PDM.xls" -r"%BaseFolder%\%%i\PDM2q_RelationsFile.xls" -f0.01 -V"PDM2q_inOuts_infoFile.txt"
set sanxotPDM2q_no_outs=sanxot.exe -aPDM2q_noOuts -p"%BaseFolder%\%%i" -d"PDM.xls" -r"%BaseFolder%\%%i\PDM2qOuts_tagged.xls" -o"protein.xls" -g -f -V"PDM2q_inOuts_infoFile.txt" --tags="!out"

REM REM REM advanced WF
REM protein site changes after protein correction
set sanxotq_PDM2qdna_in_outs=sanxot.exe -aq_PDM2qdna_inOuts -p"%BaseFolder%\%%i" -d"PDM2q_noOuts_lowerNormV.xls" -r"%BaseFolder%\%%i\PDM2qdna_RelationsFile.xls" -g --tags=""
set sanxotsieveq_PDMqdna=sanxotsieve.exe -aq_PDM2qdnaOuts -p"%BaseFolder%\%%i" -d"PDM2q_noOuts_lowerNormV.xls" -r"%BaseFolder%\%%i\PDM2qdna_RelationsFile.xls" -f0.01 -V"q_PDM2qdna_inOuts_infoFile.txt"
set sanxotq_PDM2qdna_no_outs=sanxot.exe -aq_PDM2qdna_noOuts -p"%BaseFolder%\%%i" -d"PDM2q_noOuts_lowerNormV.xls" -r"%BaseFolder%\%%i\q_PDM2qdnaOuts_tagged.xls" -o"qdna_q.xls" -g -f -V"q_PDM2qdna_inOuts_infoFile.txt" --tags="!out"

set sanxotq_qdna2qna_in_outs=sanxot.exe -aq_qdna2qna_inOuts -p"%BaseFolder%\%%i" -d"qdna_q.xls" -r"%BaseFolder%\%%i\qdna2qna_RelationsFile.xls" -g --tags=""
set sanxotsieveq_qdnaqna=sanxotsieve.exe -aq_qdna2qnaOuts -p"%BaseFolder%\%%i" -d"qdna_q.xls" -r"%BaseFolder%\%%i\qdna2qna_RelationsFile.xls" -f0.01 -V"q_qdna2qna_inOuts_infoFile.txt"
set sanxotq_qdna2qna_no_outs=sanxot.exe -aq_qdna2qna_noOuts -p"%BaseFolder%\%%i" -d"qdna_q.xls" -r"%BaseFolder%\%%i\q_qdna2qnaOuts_tagged.xls" -o"qna_q.xls" -g -f -V"q_qdna2qna_inOuts_infoFile.txt" --tags="!out"
set sanxotq_qna2A=sanxot.exe -aq_qna2A -p"%BaseFolder%\%%i" -d"qna_q.xls" -C -g --tags=""
set sanxotq_qdna2A=sanxot.exe -aq_qdna2A -p"%BaseFolder%\%%i" -d"qdna_q.xls" -C -g --tags=""

for %%i in (126 127_N 127_C 128_N 128_C 129_N 129_C 130_N 130_C 131) do (
	start /affinity 7FFF "CALCULATION_%%i" cmd.exe /K "(%aljamiaSData% & %aljamiaS2PDMRels% & %aljamiaPDM2qdnaRels% & %aljamiaqdna2qnaRels% & %aljamiaqna2qRels% & %aljamiaPDM2qRels%)"
	)
	
:wait_loop1
	for %%i in (126 127_N 127_C 128_N 128_C 129_N 129_C 130_N 130_C 131) do (
	if not exist "%BaseFolder%\%%i\PDM2q_RelationsFile.xls" goto wait_loop1
	)
  
	for %%i in (126 127_N 127_C 128_N 128_C 129_N 129_C 130_N 130_C 131) do (
	start /affinity 7FFF "CALCULATION_%%i" cmd.exe /K "(%klibrate% & %sanxotS2PDM_in_outs% & %sanxotsieveSPDM% & %sanxotS2PDM_no_outs% & %sanxotPDM2q_in_outs% & %sanxotsievePDMq% & %sanxotPDM2q_no_outs%)"
	)
	
:wait_loop2
	for %%i in (126 127_N 127_C 128_N 128_C 129_N 129_C 130_N 130_C 131) do (
	if not exist "%BaseFolder%\%%i\PDM2q_no_outs_outStats.xls" goto wait_loop2
	)
  
 	for %%i in (126 127_N 127_C 128_N 128_C 129_N 129_C 130_N 130_C 131) do (
	if not exist "%BaseFolder%\%%i" md "%BaseFolder%\%%i"
	start /affinity 7FFF "CALCULATION_%%i" cmd.exe /K "(%sanxotq_PDM2qdna_in_outs% & %sanxotsieveq_PDMqdna% & %sanxotq_PDM2qdna_no_outs% & %sanxotq_qdna2qna_in_outs% & %sanxotsieveq_qdnaqna% & %sanxotq_qdna2qna_no_outs% & %sanxotq_qna2A% & %sanxotq_qdna2A%)"
	)
	 
	
