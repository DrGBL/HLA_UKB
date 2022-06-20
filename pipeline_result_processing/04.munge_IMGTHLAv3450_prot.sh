#for this to work, need to download the alignment folder from here: https://github.com/ANHIG/IMGTHLA/tree/3450/alignments
#specifically we'll need the _prot.txt files
#you can do this by cloning the git with "git clone --branch 3450 https://github.com/ANHIG/IMGTHLA.git"
#I will now assume that you are working in the cloned git directory (which should be called "IMGTHLA")

#go to your working directory
cd /path/to/IMGTHLA/

#create a folder with temporary alignment files
mkdir -p tmp_aa_align

gene_array=("A" "B" "C" "DMA" "DMB" "DOA" "DOB" "DPA1" "DPA2" "DPB1" "DQA1" "DQB1" "DRA" "DRB" "E" "F" "G" "H" "J" "K" "L" "Y")

for gene in "${gene_array[@]}"
do
  if test -f alignments/${gene}_prot.txt; then
  
    awk '!/#/ {print}' alignments/${gene}_prot.txt | \
      tail -n +2 | \
      head -n -3 | \
      sed 's/^\s*$/new_line/g' | \
      sed 's/^\s//g' | \
      awk -v gene=${gene} 'BEGIN{RS="new_line"}{ f = "tmp_aa_align/"gene NR ".txt"; print > f; close(f) }'
      
    pos=$(head -n 1 tmp_aa_align/${gene}1.txt | \
      grep -o .*1 | \
      wc -c  | \
      awk '{print $1-1}')

    head -n 3 tmp_aa_align/${gene}1.txt | \
      tail -n 1 | \
      awk -v pos=$pos '{print substr($0,1,pos)}' > tmp_aa_align/pos_string_${gene}.txt
    
    sed -i '1i\\' tmp_aa_align/${gene}1.txt
    
    sed -i -e 1,3d tmp_aa_align/${gene}*
    
  fi
done
