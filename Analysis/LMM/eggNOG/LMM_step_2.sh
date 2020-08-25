mkdir eggNOG_explainability
mkdir GRM_files

cat ind_NEW.txt | while read -r line

do
number="$line"

gzip Kinship_mat_files/Kinship_example.grm

~/gcta_1.92.2beta_mac/bin/gcta64 --grm-gz Kinship_mat_files/Kinship_example --make-grm --out GRM_files/GRM_1

~/gcta_1.92.2beta_mac/bin/gcta64 --reml --mgrm multi_grm.txt --pheno Phen_files/phen_file__$number.txt  --reml-pred-rand --out eggNOG_explainability/GRM_reml__$number


done
exit 0
