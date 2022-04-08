library(Seurat)
print('Start adding new annotation')
print(levels(Idents(seuratObjFinal)))
seuratObjFinal     <- RenameIdents(seuratObjFinal,
                                   ## ---
                                   `12` = 'MA'  , ## 'Mast cell'
                                   ## ---
                                   `0` =  'ST', ##'Stroma',
                                   `1` = 'ST', ##'Stroma',
                                   `4` = 'ST', ##'Stroma',
                                   `5` = 'ST', ##'Stroma',
                                   `7` = 'ST', ##'Stroma',
                                   `8` = 'ST', ##'Stroma',
                                   `16` = 'ST', ##'Stroma',
                                   ## ---
                                   `10` = 'CE', ##'Ciliated epithelial'
                                   ## ---
                                   `13` = 'SE', ##'Secretory epithelial'
                                   ## ---
                                   `9` = 'P/V', ## 'Pericyte/VSMC'
                                   ## ---
                                   `2` = 'SM', ##'Smooth muscle'
                                   ## ---
                                   `6` = 'EN', ##'Endothelial'
                                   ## ---
                                   `14` = 'LE',  ##'Lymphatic endothelial
                                   ## ---
                                   `11` = 'T/NK', ## 'T/NK cell'
                                   `15` = 'T/NK', ## 'T/NK cell'
                                   `17` = 'T/NK', ## 'T/NK cell'
                                   ## ---
                                   `3` = 'MP', ## 'Macrophages'
                                   ## ---
                                   `18` = 'B/P' ## 'B/plasma cell'
)

perClusterOrder  <- c('MA', "ST", "P/V", 'SM',
                      'CE', 'SE',
                      "EN", 'LE',
                      'T/NK', 'MP', 'B/P' )

print('END adding new annotation')
