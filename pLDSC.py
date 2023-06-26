
#Add an annotation column(eSNPs or sSNPs) in baseline model
for i in {1..22}; do cat ~/database/LDSR/1000G_EUR_Phase3_baseline/"$i".txt  |awk 'BEGIN{OFS="\t"}NR==FNR&&$1=='$i'{a[NR]=$2;b[NR]=$3}NR>FNR&&FNR==1{print "scz_male_mQTLs_200kb"}NR>FNR&&FNR>=2{flag=0;for(i in a){if(a[i]<$1&&b[i]>$1){flag=1}};print flag}' scz_male_mQTLs_200kb.txt - >baselineLD."$i".txt; done

for i in {1..22}; do zcat ~/database/LDSR/1000G_EUR_Phase3_baseline/baseline."$i".annot.gz | awk 'BEGIN{OFS="\t"}NR==FNR{a[NR]=$0}NR>FNR{print a[FNR],$0}' - baselineLD."$i".txt | gzip > baselineLD."$i".annot.gz; done 

#Calculated LD score
for i in {1..22}; do ldsc.py --l2 --bfile /mnt/NAS-PD/1000Genomes/Genotypes/1000G_EUR_Phase3_plink/1000G.EUR.QC."$i" --ld-wind-cm 1 --annot baselineLD."$i".annot.gz --print-snps ~/softwares/ldsc/hapmap3_snps/hm."$i".snp --out baselineLD."$i"; done

#Make sumstats file
munge_sumstats.py --sumstats /mnt/NAS-PD/shareData/PGC/2021_PGC3_SCZ.gz --N 65967 --out PGC3.SCZ --merge- alleles ~/database/LDSR/w_hm3.snplist --a1-inc

#Partitioned heritability
ldsc.py --h2 PGC3.SCZ.sumstats.gz --ref-ld-chr baselineLD. --w-ld-chr ~/database/LDSR/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. --overlap-annot --frqfile-chr /mnt/NAS-PD/1000Genomes/Genotypes/1000G_Phase3_frq/1000G.EUR.QC. --out SCZgwas_mdmp_mQTLs_200kb --n-blocks 1000



