#1 Remove plants & unknown [using Macqiime]
filter_taxa_from_otu_table.py -i OTU_Table.Merged.99.biom -o OTU_Table.Merged.99.justFungi.biom -p k__Fungi

#2 Convert to text [using Macqiime]
biom convert -i OTU_Table.Merged.99.justFungi.biom -o OTU_Table.Merged.99.justFungi.txt -b --header-key taxonomy

#3 Eliminate bleedover (change 1's, 2's,... 9's to 0) ***Use $(printf '\t') instead of [Cont-V, then Tab] to insert tabs
cat OTU_Table.Merged.99.justFungi.txt | sed -e "s/$(printf '\t')1$(printf '\t')/$(printf '\t')0$(printf '\t')/g" | sed -e "s/$(printf '\t')2$(printf '\t')/$(printf '\t')0$(printf '\t')/g" | sed -e "s/$(printf '\t')3$(printf '\t')/$(printf '\t')0$(printf '\t')/g" | sed -e "s/$(printf '\t')4$(printf '\t')/$(printf '\t')0$(printf '\t')/g" | sed -e "s/$(printf '\t')5$(printf '\t')/$(printf '\t')0$(printf '\t')/g" | sed -e "s/$(printf '\t')6$(printf '\t')/$(printf '\t')0$(printf '\t')/g" | sed -e "s/$(printf '\t')7$(printf '\t')/$(printf '\t')0$(printf '\t')/g" | sed -e "s/$(printf '\t')8$(printf '\t')/$(printf '\t')0$(printf '\t')/g" | sed -e "s/$(printf '\t')9$(printf '\t')/$(printf '\t')0$(printf '\t')/g" > OTU_Table.Merged.99.justFungi.NoBleed.txt

#4 Convert back to biom [using Macqiime]
biom convert -i OTU_Table.Merged.99.justFungi.NoBleed.txt -o OTU_Table.Merged.99.justFungi.NoBleed.biom --table-type "otu table" --process-obs-metadata taxonomy

#5 Eliminate low read samples and Rarefy to 3000 [using Macqiime]
single_rarefaction.py -i OTU_Table.Merged.99.justFungi.NoBleed.biom -o OTU_Table.Merged.99.justFungi.NoBleed.r3K.biom -d 3000

#6 Merge species [using Macqiime]
summarize_taxa.py -i OTU_Table.Merged.99.justFungi.NoBleed.r3K.biom -o OTU_Table.M.99.jF.NB.r3K.species -L 7 -a

#7 Get out of folder
mv OTU_Table.M.99.jF.NB.r3K.species/OTU_Table.Merged.99.justFungi.NoBleed.r3K_L7.txt .

#8 Remove 0's
cat OTU_Table.Merged.99.justFungi.NoBleed.r3K_L7.txt | sed -e "s/$(printf '\t')0\.0/$(printf '\t')/g" > OTU_Table.Merged.99.justFungi.NoBleed.r3K.species.no0s.txt
