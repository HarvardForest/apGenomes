source('apg_dataloader.R')

### Compile data
temp <- read.table('../docs/Invertebrate.1.0.tsv',sep = '\t',header = T)
temp <- cbind(temp,replicate=numeric(0))
temp.cols <- colnames(temp)
temp <- do.call(rbind,lapply(1:nrow(sample.info),function(times) rep(x = '',times = length(temp.cols))))
colnames(temp) <- temp.cols
lat.lon <- apg.geo[match(sample.info[,'broadID'],rownames(apg.geo)),2:1]
orginfo <- data.frame(X.sample_name = sample.info[,'broadID'],
                      X.organism = paste('Aphaenogaster',sample.info[,'spec_epithet']),
                      X.collection_date = rep('25-Sep-2015',nrow(temp)),
                      lat_lon = paste(lat.lon[,1],lat.lon[,2],sep = '_'),
                      X.geo_loc_name = c('USA: Keystone Heights, FL',
                          'USA: Greenfield, MA',
                          'USA: Kite, GA',
                          'USA: Alachua County, FL',
                          'USA: Duke Forest, NC',
                          'USA: Duke Forest, NC',
                          'USA: Keystone Heights, FL'),
                      X.tissue = rep('whole organism',nrow(temp)),
                      breed = rep('not applicable',nrow(temp)),
                      host = rep('not applicable',nrow(temp)),
                      isolation_source = rep('not applicable',nrow(temp)),
                      isolate = rep('not applicable',nrow(temp)),
                      replicate = paste('colony',c(1,1,1,1,1,2,1),sep = '_')
)

### Rename samples
orginfo[,'X.sample_name'] <- paste(paste('Aphaenogaster',sample.info[,'spec_epithet'],sep = '_'),orginfo[,'replicate'],sep = '_')

### Adjust formatting for export
orginfo <- data.frame(orginfo,temp[,!colnames(temp) %in% colnames(orginfo)])
orginfo <- orginfo[order(orginfo[,'X.organism']),]
orginfo <- orginfo[,match(colnames(temp),colnames(orginfo))]
orginfo <- apply(orginfo,2,as.character)
### Fix columns
all(colnames(orginfo) == colnames(temp))
orginfo <- rbind(colnames(orginfo),orginfo)
orginfo <- apply(orginfo,2,gsub,pattern = 'X.',replacement = '*')

### Export
write.table(orginfo,file = '../docs/apg_orginfo.tsv',row.names = FALSE,sep = '\t',col.names = FALSE)


### Batch genome submission form
temp <- readLines('../docs/sample_batch_genome_accs_fsa.91965d47e287.tsv')
temp <- sapply(temp,strsplit,split = '\t')[[1]]
bio.sam <- read.table('../docs/biosample_info.tsv',sep = '\t')[-1,]
batch.sub <- array('',dim = c(nrow(orginfo[-1,]),length(temp)))
colnames(batch.sub) <- temp
colnames(batch.sub)
batch.sub[,'biosample_accession'] <- rep('',nrow(batch.sub))
batch.sub[,'sample_name'] <- orginfo[-1,'X.sample_name']
batch.sub[,"assembly_date"] <- rep('2016-12-16',nrow(batch.sub))
batch.sub[,"assembly_name"] <- c('AZXXQ','AZXXR','AZXXP','AZXXO','AZXXN','AJDMW','AZXXM')
batch.sub[,"assembly_method"] <- rep('ALLPATHS-LG',nrow(batch.sub))
batch.sub[,"assembly_method_version"] <- rep('not applicable',nrow(batch.sub))
### Get coverage stats
batch.sub[,"genome_coverage"] <- c(32,28,22,24,31,25,25)
batch.sub[,"sequencing_technology"] <- rep('Illumina Hi-Seq 2500',nrow(batch.sub))
batch.sub[,"reference_genome"] <- rep('not applicable',nrow(batch.sub))
batch.sub[,"update_for"] <- rep('not applicable',nrow(batch.sub))
batch.sub[,"bacteria_available_from"] <- rep('not applicable',nrow(batch.sub))
batch.sub[,"filename"] <- paste0(batch.sub[,"assembly_name"],'.tar.gz')
### Write to file
write.table(batch.sub,file = '../docs/apg_batchsub.tsv',row.names = FALSE,sep = '\t',col.names = TRUE)

