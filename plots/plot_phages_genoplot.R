# Script to plot phage genomes & to make CDS comparison. It reads all gbk files from actual dir
# and pairwise comparisons between genomes that want to be showed (in *.tab). 
# Fig. 3. 10.3390/ijms21093160

library(genoPlotR)
library(ape)
library(stringr)
library(plotfunctions)

gbks=list.files(path = ".", pattern = '*.gbk')
comps = list.files(path='.', pattern = '*.tab')

# Read dna segs from gbks
list_gbks = c()
for (i in 1:length(gbks)){
  sample=gsub('.gbk', '', gbks[i])
  assign(sample, read_dna_seg_from_file(gbks[i])) # To create the dna_segs objects
  var = as.name(sample) # To iterate trough variables
  list_gbks = c(list_gbks, var)
}
list_annot = c()

# Annotate dna_segs
for (i in 1:length(list_gbks)){
  annot_name = paste(list_gbks[i], '_annot', sep="")
  annot = eval(list_gbks[[i]])$product # Eval is used for variable conversion
  assign(annot_name, auto_annotate(eval(list_gbks[[i]]), locus_tag_pattern=NULL, names=annot))
  var = as.name(annot_name)
  list_annot = c(list_annot, var)
}
list_filter = c()
for (i in 1:length(list_annot)){
  annot_name = paste(list_annot[i], '_filter', sep="")
  index = which(eval(list_annot[[i]])$text != 'hypothetical protein') # Don't show hypothetical proteins
  assign(annot_name, eval(list_annot[[i]])[index,])
  var = as.name(annot_name)
  list_filter = c(list_filter, var)
}

# Assign colour codes according category (PHAGES) 
for (i in 1:length(list_gbks)){
  sample=gsub('.gbk', '', gbks[i])
  name=paste(sample, '_color', sep = '')
  df=eval(list_gbks[[i]])
  df$col = 'black' # By default, arrow outlines will be black
  df[df$product=='hypothetical protein', "fill"] = 'darkgrey'
  df[which(str_detect(df$product, 'ase|coll|gal', negate = FALSE)), 'fill']="#0b9663" # Other, eg. peptidase, phospoesterase, collar
  df[which(str_detect(df$product, 'DNA|nuclease', negate = FALSE)), 'fill']="#BADA55" # DNA replication
  df[which(str_detect(df$product, 'matur|Matur', negate = FALSE)), 'fill']="#7FE5F0" # DNA transcription
  df[which(str_detect(df$product, 'tail|Tail', negate = FALSE)), 'fill']="#420420" # Tail
  df[which(str_detect(df$product, 'spike', negate = FALSE)), 'fill']="purple" # Tail spike
  df[which(str_detect(df$product, 'head|int|scaf|caps|teg|Head|Int|Scaf|Cap|Teg', negate = FALSE)), 'fill']="#FF7373" # Head, capsid
  df[which(str_detect(df$product, 'lysin|spanin|holin|Lysin|Spanin|Holin', negate = FALSE)), 'fill']="#FF7F50" # Lysis
  df[which(str_detect(df$product, 'RNA', negate = FALSE)), 'fill']="#DCEDC1" # Transcription
  assign(name, df)
}

# Remove 'putative' from annotations
vB_KpnP_SU552A_annot_filter$text = gsub('putative ', '', vB_KpnP_SU552A_annot_filter$text)
VLC5_annot_filter$text = gsub('putative ', '', VLC5_annot_filter$text)
VLC6_annot_filter$text = gsub('putative ', '', VLC6_annot_filter$text)
vB_KpnP_SU552A_annot_filter$rot = 30 # We change the text angle in order to improve visualization
VLC5_annot_filter$rot = 30
VLC6_annot_filter$rot = 30
klebsiella_phage_kpv74_annot_filter$text = gsub('putative ', '', klebsiella_phage_kpv74_annot_filter$text)
klebsiella_phage_kpv74_annot_filter$rot = 30

# Read blast comparison files in proper order
comp1=try(read_comparison_from_blast("comp1_new.tab"))   # KPV74 vs VLC5                           
comp2=try(read_comparison_from_blast("comp2_new.tab"))
comp3=try(read_comparison_from_blast("comp3_new.tab")) # VLC6 vs SU5502
all_comps=rbind(comp1, comp2,comp3)

# Don't show annotation. Change if you want to show any. 
vB_KpnP_SU552A_annot_filter=NULL
klebsiella_phage_kpv74_annot_filter=NULL
VLC6_annot_filter=NULL
VLC5_annot_filter=NULL

# Plot 
plot=plot_gene_map(dna_segs = list(vB_KpnP_SU552A_color, VLC5_color, VLC6_color, klebsiella_phage_kpv74_color),
              gene_type = "arrows",dna_seg_scale = T, scale=T, annotations = list(vB_KpnP_SU552A_annot_filter, VLC5_annot_filter, VLC6_annot_filter, klebsiella_phage_kpv74_annot_filter),
              annotation_cex = 1, annotation_height=0,  
              plot_new = T, comparisons = list(comp1,comp2, comp3), legend = T, dna_seg_labels = c("SU552A","VLC5","VLC6","KpV74"),
              dna_seg_label_cex = 1.2, main = '', global_color_scheme=c("per_id", "auto", "grey", "0.5"))          

# Make %id legend
plot.new()
palette=apply_color_scheme(sort(all_comps$per_id, decreasing = F),color_scheme = "grey", transparency = 0.5)
gradientLegend(sort(all_comps$per_id, decreasing = T), color = palette, side = 1, dec = 0)

# Make CDS type & colour legend
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend('bottom', c("DNA replication","DNA packaging","Head", "Lysis", 'Other', "Tail", "Transcription", "Unknown"),
       fill=c('#BADA55', '#7FE5F0', '#FF7373', '#FF7F50', 'purple', '#420420', '#DCEDC1', 'darkgrey'))
       
