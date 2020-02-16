#!/bin/bash
#Made by Mykhaylo Usyk MSci. 2016 mu408@nyu.edu
#Annotate by Shuzhen Li. 2018.12.16

while read p
do

#splitting the unite id from the unite taxonomy file
#取出物种id。文件及最后一行的那个文件。需要提前下载。
silva_id=$(echo ${p} | cut -f 1 -d ' ')
#取出物种信息(界门纲目科属种)
silva_name=$(echo ${p} | cut -f 2 -d ' ')

echo "The unite ID was identified as ${silva_id}"

#Getting the genus name out from the unite taxonomy file
#从所有物种信息中取出属的信息
ninth_name=$(echo ${silva_name} | cut -f 6 -d ';')

#Making sure that unidentified entries won't be processed by allowing the loop to break 
#去掉未鉴定的属
if [ $ninth_name == "g__unidentified" ]
then
   echo "No identification for this genus, trying to exit" && continue
else
   echo "Proceeding..."
fi

#Removing the genus tag to make a grep id for silva
#去掉属前面的tag:g_
for_silva=$(echo ${ninth_name} | sed 's/g__//1')

#Getting the Silva id to pull out the sequence
#以上得到的UNITE中所有属信息，和silva数据库物种信息比对，挑出silva数据库中共有的属。silva_fungal_taxa.txt需提前下载。
silva_entry=$(grep ";${for_silva}[$;]" ./silva_fungal_taxa.txt)

if [ -z "$silva_entry" ]
then
   echo "A matching Silva entry to ${for_silva} wasn't found, moving on" && continue
fi

echo "Matching Silva entry was found, trying to concatenate 18S with unite ITS entry"

silva_ref=$(echo ${silva_entry} | head -n 1 | cut -f 1 -d ' ')

#Getting the silva sequence by using the obtained id
#从silva数据库中挑出共有属的序列。文件需下载。
eis_seq=$(awk "/>${silva_ref}$/{getline; print}" ./silva_fungal_joined.fasta)
#从UNITE数据库中挑出共有属的序列。文件需下载。
ITS_seq=$(awk "/${silva_id}/{getline; print}" ./sh_refs_qiime_ver7_97_31.01.2016.fasta)

#Time to write the sequence to file, but a final check of all the variable so as not to mess up the output
if [ -z "$eis_seq" ]
then
   echo "NO 18S SEQUENCE, STOPPING" && continue
fi

if [ -z "$ITS_seq" ]
then
   echo "NO ITS SEQUENCE, STOPPING" && continue
fi

#输出。第一行为UNITE数据库中的序列id，及属水平注释信息。第二行为18S序列+ITS序列。
echo ">${silva_id}" >> great_pizza_toppings.txt
echo "${eis_seq}${ITS_seq}" >> great_pizza_toppings.txt


done < sh_taxonomy_qiime_ver7_97_31.01.2016.txt
