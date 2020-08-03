mkdir Kegg_explainability
#mkdir Kegg_explainability_FF
mkdir GRM_files

cat ind_NEW.txt | while read -r line

do
number="$line"

gzip Kinship_mat_files/Kinship_env.grm
#gzip Kinship_mat_files/Kinship_genes.grm

~/gcta_1.92.2beta_mac/bin/gcta64 --grm-gz Kinship_mat_files/Kinship_env --make-grm --out GRM_files/GRM_1
#~/gcta_1.92.2beta_mac/bin/gcta64 --grm-gz Kinship_mat_files/Kinship_genes --make-grm --out GRM_files/GRM_2

#~/gcta_1.92.2beta_mac/bin/gcta64 --reml --mgrm multi_grm.txt --pheno Phen_files/phen_file__1802.txt  --reml-pred-rand --qcovar Fixed_effects/Fixed_effects.txt --reml-est-fix --out Kegg_explainability_FF/GRM_reml__1

~/gcta_1.92.2beta_mac/bin/gcta64 --reml --mgrm multi_grm.txt --pheno Phen_files/phen_file__$number.txt  --reml-pred-rand --out Kegg_explainability/GRM_reml__$number


done
exit 0


#~/gcta_1.92.2beta_mac/bin/gcta64 --grm-gz Kinship_mat_files/Kinship_exampl --make-grm --out GRM_files/GRM_1
#
#~/gcta_1.92.2beta_mac/bin/gcta64 --reml --mgrm multi_grm.txt --pheno Phen_files/phen_file__476.txt  --reml-pred-rand --qcovar Fixed_effects/Fixed_effects.txt --reml-est-fix --out Kegg_explainability_FF/GRM_reml__1
#
#
#~/gcta_1.92.2beta_mac/bin/gcta64 --reml --mgrm multi_grm.txt --pheno Phen_files/phen_file__476.txt  --reml-pred-rand --out Kegg_explainability/GRM_reml__1
