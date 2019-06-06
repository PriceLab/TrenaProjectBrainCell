library(TrenaProjectBrainCell)
library(RUnit)
library(trenaSGM)
library(org.Hs.eg.db)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tp")) {
   message(sprintf("--- creating instance of TrenaProjectBrainCell"))
   tp <- TrenaProjectBrainCell();
}
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()
   test_supportedGenes()
   test_variants()
   test_footprintDatabases()
   test_expressionMatrices()
   test_setTargetGene()
   test_buildSingleGeneModel()
   test_buildSingleGeneModel_slowGenes()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   message(sprintf("--- test_constructor"))

   checkTrue(all(c("TrenaProjectBrainCell", "TrenaProjectHG38") %in% is(tp)))
   checkEquals(getFootprintDatabasePort(tp), 5432)

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_supportedGenes <- function()
{
   message(sprintf("--- test_supportedGenes"))

   subset.expected <- c("APOE")
   checkTrue(all(subset.expected %in% getSupportedGenes(tp)))

} # test_supportedGenes
#------------------------------------------------------------------------------------------------------------------------
test_variants <- function()
{
   message(sprintf("--- test_variants"))

   checkEquals(length(getVariantDatasetNames(tp)), 0)

} # test_variants
#------------------------------------------------------------------------------------------------------------------------
test_footprintDatabases <- function()
{
   message(sprintf("--- test_footprintDatabases"))

   expected <- c("brain_hint_16", "brain_hint_20", "brain_wellington_16", "brain_wellington_20")

   checkTrue(all(expected %in% getFootprintDatabaseNames(tp)))
   checkEquals(getFootprintDatabaseHost(tp), "khaleesi.systemsbiology.net")
   checkEquals(getFootprintDatabasePort(tp), 5432)

} # test_footprintDatabases
#------------------------------------------------------------------------------------------------------------------------
test_expressionMatrices <- function()
{
   message(sprintf("--- test_expressionMatrices"))
    
   expected.matrix <- "Micro_TYROBP"
   checkTrue(all(expected.matrix %in% getExpressionMatrixNames(tp)))

   mtx <- getExpressionMatrix(tp, expected.matrix)
   checkEquals(dim(mtx), c(50281,  264))
   expected.genes <- c("A1BG", "A1CF", "A2M")
   checkTrue(all(expected.genes %in% rownames(mtx)))

   checkTrue(max(mtx) < 15)

} # test_expressionMatrices
#------------------------------------------------------------------------------------------------------------------------
# setting the target gene implies a few other assignements, all tested here:
#   geneInfo (temporarily also masquerading at tbl.transcripts
#   geneRegion
#   geneEnhancersRegion (when avaialable, defaults to geneRegion)
#
test_setTargetGene <- function()
{
   message(sprintf("--- test_setTargetGene"))

   setTargetGene(tp, "MICA")
   checkEquals(getTargetGene(tp), "MICA")

   message(sprintf("    transcripts"))
   tbl.transcripts <- getTranscriptsTable(tp)
   checkTrue(nrow(tbl.transcripts) == 1)
   checkEquals(tbl.transcripts$chr, "chr6")

   checkEquals(tbl.transcripts$start, 31399784)
   checkEquals(tbl.transcripts$end , 31415315)
   checkEquals(tbl.transcripts$tss, 31403579)
   checkEquals(tbl.transcripts$strand, 1)

   message(sprintf("    geneRegion"))
   region <- getGeneRegion(tp, flankingPercent=0)
   checkTrue(all(c("chromLocString", "chrom", "start", "end") %in% names(region)))
   checkEquals(region$chromLocString, "chr6:31399784-31415315")

   message(sprintf("    enhancers"))
   tbl.enhancers <- getEnhancers(tp)
   checkEquals(colnames(tbl.enhancers), c("chrom", "start", "end", "type", "combinedScore", "geneSymbol"))
   checkTrue(nrow(tbl.enhancers) >= 0)

   message(sprintf("    geneGeneEnhancersRegion"))
   region <- getGeneEnhancersRegion(tp, flankingPercent=0)
   checkTrue(all(c("chromLocString", "chrom", "start", "end") %in% names(region)))
   checkEquals(region$chromLocString, "chr6:30554823-32372101")

   message(sprintf("    encode DHS"))
   tbl.dhs <- getEncodeDHS(tp)
   checkTrue(nrow(tbl.dhs) > 1900)

   message(sprintf("    ChIP-seq"))
   tbl.chipSeq <- with(tbl.transcripts, getChipSeq(tp, chrom=chrom, start=start, end=end, tfs="BCLAF1"))
   checkEquals(nrow(tbl.chipSeq), 2)

} # test_setTargetGene
#------------------------------------------------------------------------------------------------------------------------
test_buildSingleGeneModel <- function()
{
   printf("--- test_buildSingleGeneModel")

   genome <- "hg38"
   targetGene <- "APOE"
   chromosome <- "chr19"
   tss <- 44905751
      # strand-aware start and end: trem2 is on the minus strand
   geneHancer.promoter.chromLocString <- "chr19:44,903,353-44,907,298"
   start <- 44903353
   end   <- 44907298
   tbl.regions <- data.frame(chrom=chromosome, start=start, end=end, stringsAsFactors=FALSE)
   matrix.name <- "Micro_TYROBP"
   checkTrue(matrix.name %in% getExpressionMatrixNames(tp))
   mtx <- getExpressionMatrix(tp, matrix.name)

   build.spec <- list(title="unit test on APOE",
                      type="footprint.database",
                      regions=tbl.regions,
                      geneSymbol=targetGene,
                      tss=tss,
                      matrix=mtx,
                      db.host=getFootprintDatabaseHost(tp),
                      db.port=getFootprintDatabasePort(tp),
                      databases=getFootprintDatabaseNames(tp),
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      motifDiscovery="builtinFimo",
                      tfPool=allKnownTFs(),
                      tfMapping="MotifDB",
                      tfPrefilterCorrelation=0.1,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   fpBuilder <- FootprintDatabaseModelBuilder(genome, targetGene,  build.spec, quiet=FALSE)
   suppressWarnings(x <- build(fpBuilder))
   lapply(x, dim)
   
   checkEquals(sort(names(x)), c("model", "regulatoryRegions"))
   tbl.regulatoryRegions <- x$regulatoryRegions
   tbl.model <- x$model
   tbl.model <- tbl.model[order(tbl.model$rfScore, decreasing=TRUE),]
   checkTrue(all(tbl.model$gene %in% tbl.regulatoryRegions$geneSymbol))
   checkTrue(nrow(x$model) > 50)
   checkTrue("TP53" %in% head(x$model$gene))
   checkTrue(max(tbl.model$pearsonCoeff) > 0.85)
     # a modest sanity check on pearsonCoeff: should be exactly what we see in the expression matrix
   checkEqualsNumeric(cor(mtx["APOE",], mtx["TP53",]), subset(tbl.model, gene=="TP53")$pearsonCoeff)

} # test_buildSingleGeneModel
#------------------------------------------------------------------------------------------------------------------------
# cory's top genes producing huge output (5 jun 2019): "SF3A2" "ZNF764" "PRR12" "ALDH16A1" "EIF1AD" "ZNF44"   
test_buildSingleGeneModel_slowGenes <- function()
{
   printf("--- test_buildSingleGeneModel_slowGenes")

   slowGenes <- c("SF3A2", "ZNF764", "PRR12", "ALDH16A1", "EIF1AD", "ZNF44")
   genome <- "hg38"
   targetGene <- slowGenes[1]
   setTargetGene(tp, slowGenes[1])

   tss <- getTranscriptsTable(tp, targetGene)$tss
   
   tbl.regions <- getEnhancers(tp)[, c("chrom", "start", "end")]
   tbl.regions <- tbl.regions[order(tbl.regions$start, decreasing=FALSE),]
   matrix.name <- "Micro_TYROBP"
   checkTrue(matrix.name %in% getExpressionMatrixNames(tp))
   mtx <- getExpressionMatrix(tp, matrix.name)

   build.spec <- list(title=sprintf("unit test on %s", targetGene), 
                      type="footprint.database",
                      regions=tbl.regions,
                      geneSymbol=targetGene,
                      tss=tss,
                      matrix=mtx,
                      db.host=getFootprintDatabaseHost(tp),
                      db.port=getFootprintDatabasePort(tp),
                      databases=getFootprintDatabaseNames(tp),
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      motifDiscovery="builtinFimo",
                      tfPool=allKnownTFs(),
                      tfMapping=c("MotifDB"), # , "TFClass"),
                      tfPrefilterCorrelation=0.5,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   fpBuilder <- FootprintDatabaseModelBuilder(genome, targetGene,  build.spec, quiet=FALSE)
   print(system.time(suppressWarnings(x <- build(fpBuilder))))

   checkEquals(sort(names(x)), c("model", "regulatoryRegions"))
   tbl.regulatoryRegions <- x$regulatoryRegions
   tbl.model <- x$model
   tbl.model <- tbl.model[order(abs(tbl.model$pearsonCoeff), decreasing=TRUE),]
   checkTrue(all(tbl.model$gene %in% tbl.regulatoryRegions$geneSymbol))
   checkTrue(nrow(tbl.model) > 50)
   sample.strong.tf <- tbl.model$gene[1]
   checkTrue(max(tbl.model$pearsonCoeff) > 0.5)
     # a modest sanity check on pearsonCoeff: should be exactly what we see in the expression matrix
   checkEqualsNumeric(cor(mtx[targetGene,], mtx[sample.strong.tf,]), subset(tbl.model, gene==sample.strong.tf)$pearsonCoeff)

} # test_buildSingleGeneModel_slowGenes
#------------------------------------------------------------------------------------------------------------------------
# no genehancer info for this gene, and though we have gene expression, there are no footprints.
# do we handle this oddity gracefully?
test_buildSingleGeneModel_RBMXP2 <- function()
{
   printf("--- test_buildSingleGeneModel_RBMXP2")

   genome <- "hg38"
   targetGene <- "RBMXP2"

   chromosome <- "chr9"
   tss <- 30689105
   start <- tss - 5000
   end <- tss + 5000

   tbl.regions <- data.frame(chrom=chromosome, start=start, end=end, stringsAsFactors=FALSE)
   matrix.name <- "GTEx.liver.geneSymbols.matrix.asinh"
   checkTrue(matrix.name %in% getExpressionMatrixNames(tp))
   mtx <- getExpressionMatrix(tp, matrix.name)

   build.spec <- list(title="unit test on RBMXP2",
                      type="footprint.database",
                      regions=tbl.regions,
                      geneSymbol=targetGene,
                      tss=tss,
                      matrix=mtx,
                      db.host=getFootprintDatabaseHost(tp),
                      db.port=getFootprintDatabasePort(tp),
                      databases=getFootprintDatabaseNames(tp),
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      motifDiscovery="builtinFimo",
                      tfPool=allKnownTFs(),
                      tfMapping="MotifDB",
                      tfPrefilterCorrelation=0.1,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   fpBuilder <- FootprintDatabaseModelBuilder(genome, targetGene,  build.spec, quiet=FALSE)
   checkException(x <- build(fpBuilder), silent=TRUE)

} # test_buildSingleGeneModel_RBMXP2
#------------------------------------------------------------------------------------------------------------------------

if(!interactive())
   runTests()
