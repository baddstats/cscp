## note path relative to this file's location
load("../data/cancer_example.RData")
png("../plots/cancer_pp.png")
plot(stroma_cells$inner_x, stroma_cells$inner_y, pch = 20, axes = FALSE, xlab = "", ylab = "")
points(tumour_cells$inner_x, tumour_cells$inner_y, pch = 20, col = "red")
box()
dev.off()
