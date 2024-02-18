#filter vcf
vcftools --vcf raw.vcf \
#--remove lowDP.indv \
--max-alleles 2 \
--min-alleles 2 \
--minDP 4 \
--minQ 30 \
--mac 4  \
--max-missing 0.85 \
--out raw.2 \
--recode \
--recode-INFO-all

#annotate vcf
#java -jar spEff/snpEff.jar Hb raw.2.recode.vcf > raw.annotated.vcf

#LD filter
./plink --vcf ../raw.2.recode.vcf --maf 0.05 --indep 50 5 2 --allow-extra-chr
./plink --vcf ../raw.2.recode.vcf --extract plink.prune.in --out LD-filter --recode vcf-iid --allow-extra-chr --keep-allele-order

#Extract SNPs of Single Copy Genes 
makeblastdb -in gene.fasta -dbtype nucl -out ref
blastn -query gene.fasta -out gene.blast -db ref -outfmt 6 -evalue 1e-5 -num_threads 10
perl get_SCbed.pl gene.blast ref.gff3 out.bed
perl scVCF.pl -v raw.2.recode.vcf -b out.bed -o clean.vcf
perl VCF2fasta.v2.pl -v clean.vcf -o snp.fasta
#construct ML tree
modeltest-ng-static -i snp.fasta -d nt
raxmlHPC-PTHREADS -T 12 -s snp.fasta -m GTRCAT -n tree2 -p 12345 -f a -# 100 -x 12345
raxmlHPC-PTHREADS-SSE3 -T 50 -s snp.fasta -m GTRGAMMAX -n tree-new -p 12345 -f a -# 100 -x 12345

#PCA
sed -i 's/chr//g' raw.2.recode.vcf                                       
vcftools --vcf raw.2.recode.vcf --plink --out tmp			
./plink --noweb --file tmp --make-bed --out tmp --chr-set 7		
./gcta64 --bfile tmp --make-grm --autosome --out tmp --autosome-num 7	
./gcta64 --grm tmp --pca 3 --out pcatmp                                  

#Structure
for K in {1..9} ; do ./admixture --cv tmp.bed $K | tee log${K}.out; done

#tajimaD
vcftools --vcf treat.vcf --TajimaD 50000 --keep sample.list --out tjmD_50k
#pi
vcftools --vcf treat.vcf --keep population1.txt --window-pi 50000 --out p1-50k
vcftools --vcf treat.vcf --keep population2.txt --window-pi 50000 --out p2-50k
#CLR
for i in {1..7}; do vcftools --vcf raw.2.recode.vcf --recode --recode-INFO-all --stdout --chr ${i} --keep I.txt > I.chr${i}.vcf; done
ls *.txt | while read id; do \
~/opt/biosoft/sweed-master/SweeD-P -name $id.chr1.50kb -input $id.chr1.vcf -grid 7822 -minsnps 200 -maf 0.05 -missing 0.1 -threads 3
~/opt/biosoft/sweed-master/SweeD-P -name $id.chr2.50kb -input $id.chr2.vcf -grid 10806 -minsnps 200 -maf 0.05 -missing 0.1 -threads 3
~/opt/biosoft/sweed-master/SweeD-P -name $id.chr3.50kb -input $id.chr3.vcf -grid 10436 -minsnps 200 -maf 0.05 -missing 0.1 -threads 3
~/opt/biosoft/sweed-master/SweeD-P -name $id.chr4.50kb -input $id.chr4.vcf -grid 9522 -minsnps 200 -maf 0.05 -missing 0.1 -threads 3
~/opt/biosoft/sweed-master/SweeD-P -name $id.chr5.50kb -input $id.chr5.vcf -grid 9474 -minsnps 200 -maf 0.05 -missing 0.1 -threads 3
~/opt/biosoft/sweed-master/SweeD-P -name $id.chr6.50kb -input $id.chr6.vcf -grid 8992 -minsnps 200 -maf 0.05 -missing 0.1 -threads 3
~/opt/biosoft/sweed-master/SweeD-P -name $id.chr7.50kb -input $id.chr7.vcf -grid 10682 -minsnps 200 -maf 0.05 -missing 0.1 -threads 3;
done

#treemix
vcftools --vcf test.vcf --plink-tped --out test
plink --bfile tmp --recode --transpose --out test
plink --tfile test --freq  --within pop.cov
gzip plink.frq.strat
python plink2treemix.py plink.frq.strat.gz treemix_in.gz
treemix -se -bootstrap -k 1000 -m 1 -i treemix_in.gz
for m in {1..10};do treemix -se -bootstrap -k 1000 -m ${m} -i treemix_in.gz -global -o TreeMix.${m}; done
#ABBA
Dsuite Dtrios YourDataName-id-maf0.05.vcf sample-group.txt

#PopCNV
popCNV -g /public/home/Hordeum_genome/Hb.fasta -s 1000 -t 6 -r /public/home/Hordeum_genome/Resequencing_reads/CNV/I/read_depth/ -b /public/home/Hordeum_genome/Resequencing_reads/CNV/I/ -l /public/home/Hordeum_genome/Resequencing_reads/CNV/gene.list -w  /public/home/Hordeum_genome/Resequencing_reads/CNV/I/ --group group.list --sample sample.list --wild I1