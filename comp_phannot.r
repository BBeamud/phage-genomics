library(ape)
library(compare)
library(data.table)
library(VennDiagram)
setwd("~/Documentos/proyectos/fpu_resistencias/secuencias/phages/from_fastq_Justine/3.4_Annot_comp/")
comp=read.csv("4_1_comp.csv", header = T, dec=',')


prokka = read.gff("4-1.gff3")
prokka = prokka[,c('start','end')]
phannot = read.table("4-1_sc1_reordered_phanotate.out", col.names = c("start", "end", "frame", "contig", "score"))
phannot = phannot[,c('start', 'end', 'score')]
phaster = read.table("phaster_4.1", col.names = c("start", "end"))

# We assume PHANNOTATE is the method that detect more CDS
n_prokka = nrow(prokka)
n_phannot = nrow(phannot)
n_phaster = nrow(phaster)
a=c(n_prokka, n_phannot, n_phaster)

number= which(a==max(a))
if (number != 2){
  print("PHAHHONATE IS NOT THE METHOD WITH MORE CDS DETECTED")
  exit()
} 

prokka$id = rep(1:n_prokka)
phannot$id = rep(1:n_phannot)
phaster$id = rep(1:n_phaster)

# Function to determine the closest number in a list to other number
closest = function(number, list){
  dt = data.table(list, val=list)
  setattr(dt, "sorted", "value")
  setkey(dt, val)
  dt[J(number), roll = "nearest"][1][[1]]
}


comparison <- data.frame(PROKKA_ID=integer(), PROKKA_START=integer(), PROKKA_END=integer(), PHANNOTATE_ID=integer(),
                         PHANNOTATE_START = integer(), PHANNOTATE_END=integer(), PHANNOTATE_SCORE =integer(), 
                         PHASTER_ID = integer(), PHASTER_START=integer(), PHASTER_END=integer()) 

for (i in 1:nrow(phannot)){
  start = phannot$start[i]
  end = phannot$end[i]
  s_dif_prokka = abs(start - closest(start, prokka$start))
  e_dif_prokka = abs(end - closest(end, prokka$end))
  s_dif_phaster = abs(start - closest(start, phaster$start))
  e_dif_phaster = abs(end - closest(end, phaster$end))
  dif_prokka = s_dif_prokka + e_dif_prokka
  dif_phaster = s_dif_phaster + e_dif_phaster
  comparison[i, "PHANNOTATE_ID"] = phannot$id[i]
  comparison[i, "PHANNOTATE_START"] = phannot$start[i]
  comparison[i, "PHANNOTATE_END"] = phannot$end[i]
  comparison[i, "PHANNOTATE_SCORE"] = phannot$score[i]
  if (dif_prokka > 75) {
    comparison[i, "PROKKA_ID"] = NA
    comparison[i, "PROKKA_START"] = NA
    comparison[i, "PROKKA_END"] = NA
  } else {
    index_s = which(prokka$start %in% closest(start, prokka$start))
    index_e = which(prokka$end %in% closest(end, prokka$end))  
    if (index_s != index_e){
      print("Coordinate failled")
      exit()
    }
    comparison[i, "PROKKA_START"] = prokka$start[index_s]
    comparison[i, "PROKKA_END"] = prokka$end[index_e]
    comparison[i, "PROKKA_ID"] = prokka$id[index_s]
  }
  if (dif_phaster > 75) {
    comparison[i, "PHASTER_ID"] = NA
    comparison[i, "PHASTER_START"] = NA
    comparison[i, "PHASTER_END"] = NA
  } else {
    index_s = which(phaster$start %in% closest(start, phaster$start))
    index_e = which(phaster$end %in% closest(end, phaster$end))  
    if (index_s != index_e){
      print("Coordinate failled")
      exit()
    }
    comparison[i, "PHASTER_START"] = phaster$start[index_s]
    comparison[i, "PHASTER_END"] = phaster$end[index_e]
    comparison[i, "PHASTER_ID"] = phaster$id[index_s]
  }
}

# We need to iterate in prokka and phaster missing cds
prokka_done = unique(comparison$PROKKA_ID)
prokka_miss = which(prokka$id %in% prokka_done == FALSE)


counter=nrow(comparison)
for(i in prokka_miss){
  start = prokka$start[i]
  end = prokka$end[i]
  s_dif_phaster = abs(start - closest(start, phaster$start))
  e_dif_phaster = abs(end - closest(end, phaster$end))
  dif_phaster = s_dif_phaster + e_dif_phaster
  comparison[counter+1, "PROKKA_ID"] = prokka$id[i]
  comparison[counter+1, "PROKKA_START"] = prokka$start[i]
  comparison[counter+1, "PROKKA_END"] = prokka$end[i]
  if (dif_phaster > 75) {
    comparison[counter+1, "PHASTER_ID"] = NA
    comparison[counter+1, "PHASTER_START"] = NA
    comparison[counter+1, "PHASTER_END"] = NA
  } else {
    index_s = which(phaster$start %in% closest(start, phaster$start))
    index_e = which(phaster$end %in% closest(end, phaster$end))  
    if (index_s != index_e){
      print("Coordinate failled")
      exit()
    }
    comparison[counter+1, "PHASTER_START"] = phaster$start[index_s]
    comparison[counter+1, "PHASTER_END"] = phaster$end[index_e]
    comparison[counter+1, "PHASTER_ID"] = phaster$id[index_s]
  }
  counter = counter +1
}

phaster_done = unique(comparison$PHASTER_ID)
phaster_miss = which(phaster$id %in% phaster_done == FALSE)

for(i in phaster_miss){
  start = phaster$start[i]
  end = phaster$end[i]
  s_dif_prokka = abs(start - closest(start, prokka$start))
  e_dif_prokka = abs(end - closest(end, prokka$end))
  dif_prokka = s_dif_prokka + e_dif_prokka
  comparison[counter+1, "PHASTER_ID"] = phaster$id[i]
  comparison[counter+1, "PHASTER_START"] = phaster$start[i]
  comparison[counter+1, "PHASTER_END"] = phaster$end[i]
  if (dif_prokka > 75) {
    comparison[counter+1, "PROKKA_ID"] = NA
    comparison[counter+1, "PROKKA_START"] = NA
    comparison[counter+1, "PROKKA_END"] = NA
  } else {
    index_s = which(prokka$start %in% closest(start, prokka$start))
    index_e = which(prokka$end %in% closest(end, prokka$end))  
    if (index_s != index_e){
      print("Coordinate failled")
      exit()
    }
    comparison[counter+1, "PROKKA_START"] = prokka$start[index_s]
    comparison[counter+1, "PROKKA_END"] = prokka$end[index_e]
    comparison[counter+1, "PROKKA_ID"] = prokka$id[index_s]
  }
  counter = counter +1
}


comp=comparison

# Venn diagram
n13=nrow(subset(comp, PROKKA_ID > 0 & PHANNOTATE_ID > 0)) # PROKKA & PHANNOTATE 
n23=nrow(subset(comp, PROKKA_ID > 0 & PHASTER_ID > 0)) # PROKKA & PHASTHER
n12=nrow(subset(comp, PHANNOTATE_ID > 0 & PHASTER_ID > 0)) # PHASTER & PHANNOTATE
n123=nrow(subset(comp, PHANNOTATE_ID > 0 & PHASTER_ID > 0 & PROKKA_ID > 0)) # THE THREE
n3=nrow(subset(comp, PROKKA_ID > 0 & is.na(PHANNOTATE_ID) & is.na(PHASTER_ID))) # UNIQ PROKKA
n1=nrow(subset(comp, is.na(PROKKA_ID) & PHANNOTATE_ID > 0 & is.na(PHASTER_ID))) # UNIQ PHANNOTATE
n2=nrow(subset(comp, is.na(PROKKA_ID) & PHANNOTATE_ID > 0 & PHASTER_ID>0)) # UNIQ PHASTHER

#grid.newpage()
#draw.triple.venn(area1 = max(phannot$id), area2 = max(phaster$id), area3 = max(prokka$id), n12 = n12, n23 = n23, n13 = n13, 
#                 n123 = n123, category = c("PHANNOTATE", "GLIMMER", "PRODIGAL"), lty = "blank", 
#                 fill = c("skyblue", "pink1", "mediumorchid"))




for (i in 1:nrow(comp)) {
  prokka = (is.na(comp$PROKKA_ID[i]))
  phannot = (is.na(comp$PHANNOTATE_ID[i]))
  phaster = (is.na(comp$PHASTER_ID[i]))
  if (prokka == T){
    prokka = 0
  } else {
    prokka = 1
  }
  if (phannot == T) {
    phannot = 0
  } else {
    phannot = 1
  }
    if (phaster == T){
      phaster = 0
    } else {
      phaster = 1
    } 
  punct = phannot + prokka + phaster
  # We also considered phannotate events detected only with phannotate if score <= -3
  x = comp$PHANNOTATE_SCORE[i]
  if ((!is.na(x)) & x <= -3){
    punct = punct + 0.5
  }
  comp$punct[i] = punct
}

three = which(comp$punct >= 3)
two =  which(comp$punct >= 2 & comp$punct < 3)
one = which(comp$punct >= 1.5 & comp$punct < 2)
comp$start = 0
comp$end = 0


for(i in three){
  # We trust prokka and phaster vs phannotate 
  equal = length(unique(c(comp$PHASTER_START[i], comp$PROKKA_START.[i]))) == 1 
  if (equal == T){
    comp$start[i] = comp$PHASTER_START[i]
  } else {
    
  }
  equal = length(unique(c(comp$PHASTER_END[i], comp$PROKKA_END[i]))) == 1 
  if (equal == T){
    comp$end[i] = comp$PHASTER_END[i]
  } else {
    
  }
}

for(i in two){
  # First we need to determine which methods detected 
  prokka = is.na(comp$PROKKA_ID[i])
  phaster = is.na(comp$PHASTER_ID[i])
  phannot = is.na(comp$PHANNOTATE_ID[i])
  if (prokka != TRUE){
    comp$start[i] = comp$PROKKA_START[i]
    comp$end[i] = comp$PROKKA_END[i]
  } 
  if (phaster != TRUE) {
    comp$start[i] = comp$PHASTER_START[i]
    comp$end[i] = comp$PHASTER_END[i]
  }
}

for(i in one){
  comp$start[i] = comp$PHANNOTATE_START[i]
  comp$end[i] = comp$PHANNOTATE_END[i]
}

summary=comp[comp$punct>=1.5,c("start", "end")]


summary_sort = summary[order(summary$start, summary$end),]
write.table(summary_sort, "summary_sort", row.names = F)
