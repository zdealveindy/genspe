#' Tcltk interfaces for \code{calculate.theta} function
#' @author David Zeleny
#' @param input.matrix Matrix (or data.frame) with species composition data.
#' @param species.data Matrix (or data.frame) with species data. If suplied, it should have at least two columns - the first containing species name, the second containing layer. 
#' @param juicer Argument specific for launching the function from JUICE software; logical (default = F) - is the function launched from JUICE? If \code{juicer = TRUE}, function is expecting that \code{species.data} have JUICE-specific structure, which enables to import data back to JUICE.

#' @details Function \code{calculate.theta.tcltk} launches tcltk clickable interface, which enables to select methods and parameters used for calculation. 

#' @export
calculate.theta.tcltk <- function (input.matrix, species.data = NULL, juicer = T)
{
require (tcltk)
cancel <- tclVar (0)
end.end <- F
beals.file <- NULL

GSmethods <- c ("Additive beta diversity (Fridley et al. 2007)", "Multiplicative beta diversity (Zeleny 2009)", "Multiplicative beta on species pool (Botta-Dukat 2012)", "Pairwise Jaccard dissimilarity (Manthey & Fridley 2009)", "Multiple Sorensen dissimilarity (Manthey & Fridley (2009)")
base <- tktoplevel()
tkwm.title(base, "Generalists-specialists")

spec.frm <- tkframe (base, borderwidth=2)
frame.title <- tkframe (spec.frm, relief = 'groove', borderwidth = 2, background = 'grey')
frame.a <- tkframe (spec.frm, relief = 'groove', borderwidth = 2)
frame.b <- tkframe (spec.frm, relief = 'groove', borderwidth = 2)
frame.c <- tkframe (spec.frm, relief = 'groove', borderwidth = 2)
frame.d <- tkframe (spec.frm, borderwidth = 2)
frame.e <- tkframe (spec.frm, relief = 'groove', borderwidth = 2)
frame.e.1 <- tkframe (frame.e)
frame.e.2 <- tkframe (frame.e)
frame.f <- tkframe (spec.frm, relief = 'groove', borderwidth = 2)
frame.f.1 <- tkframe (frame.f)
frame.g <- tkframe (spec.frm, relief = 'groove', borderwidth = 2)


GSmethod <- tclVar ("additive")
label.radio <- tklabel (frame.a, text = "Which beta diversity algorithm to use?")
radio1 <- tkradiobutton (frame.a, text = GSmethods[1], value = "additive", variable = GSmethod)
radio2 <- tkradiobutton (frame.a, text = GSmethods[2], value = "multiplicative", variable = GSmethod)
radio3 <- tkradiobutton (frame.a, text = GSmethods[3], value = "beals", variable = GSmethod)
radio4 <- tkradiobutton (frame.a, text = GSmethods[4], value = "pairwise.jaccard", variable = GSmethod)
radio5 <- tkradiobutton (frame.a, text = GSmethods[5], value = "multi.sorensen", variable = GSmethod)
tk.psample <- tclVar (5)
tk.reps <- tclVar (10)
parallel <- tclVar (0)
no.cores <- tclVar (2)
remove.out <- tclVar (0)
label.entry1 <- tklabel (frame.b, text = "Minimal frequency of species ")
entry1 <- tkentry (frame.b, width = 5, textvariable = tk.psample)

label.entry2 <- tklabel (frame.c, text = "Number of random subsamples ")
entry2 <- tkentry (frame.c, width = 5, textvariable = tk.reps)

button1 <- tkbutton (frame.d, text = "Calculate", width = 10, height = 2, command = function () end.end <- calculate.theta (input.matrix = input.matrix, species.data = species.data, psample = as.numeric (tclvalue (tk.psample)), reps = as.numeric (tkget (entry2)), method = as.character (tclvalue (GSmethod)), beals.file = beals.file, parallel = as.logical (as.numeric (tclvalue (parallel))), no.cores = as.numeric (tclvalue (no.cores)), remove.out = as.logical (as.numeric (tclvalue (remove.out))), verbal = T, juicer = T, tcltk = T))


choose.label <- tklabel (frame.e.2, text = 'Select the file with beals smoothed data')
choose.button <- tkbutton (frame.e.1, text = 'Select', command = function () assign ('beals.file', choose.files (), env = .GlobalEnv))
tkpack (choose.button)
tkpack (choose.label)
tkpack (tklabel (frame.e, text = 'Beals smoothing (included in method of Botta-Dukat 2012)'), anchor = 'w')
tkpack (frame.e.1, frame.e.2, side = 'left',ipady = 5, ipadx = 5, padx = 5, pady = 5)


parallel.label <- tklabel (frame.f, text = 'Parallel calculation (enable only if you have more than one core)')
parallel.no.cores.label <- tklabel (frame.f.1, text = 'number of cores: ')
parallel.no.cores.entry <- tkentry (frame.f.1, width = 2, textvariable = no.cores)
parallel.checkbutton <- tkcheckbutton (frame.f.1, text = 'use parallel computing,', variable = parallel)

tkpack (tklabel (frame.g, text = 'Outlier analysis (McCune & Mefford 1999, suggested by Botta-Dukat 2012)'), tkcheckbutton (frame.g, text = 'remove outlier samples (with very different species composition)', variable = remove.out), anchor = 'w')

tkpack (label.radio, radio1, radio2, radio4, radio5, radio3, anchor = 'w')
tkpack (label.entry1, entry1, anchor = 'w', side = 'left')
tkpack (label.entry2, entry2, anchor = 'w', side = 'left')
tkpack (button1)
tkpack (parallel.checkbutton, parallel.no.cores.label, parallel.no.cores.entry, side = 'left')
tkpack (parallel.label,  frame.f.1, anchor = 'w')

tkpack (tklabel (frame.title, text = 'Calculation of generalists and specialists using co-occurrence species data \n Author: David Zeleny (zeleny.david@gmail.com) \n JUICE-R application (www.bit.ly/habitat-specialists)'), ipady = 10, ipadx = 10, padx = 10, pady = 10)

tkpack (frame.title, side = 'top', expand = T, fill = 'both')
tkpack (frame.a, side = "top", ipady = 10, ipadx = 10, padx = 10, pady = 10, anchor = "w", expand = T, fill = 'both')
tkpack (frame.e, ipady = 10, ipadx = 10, padx = 10, pady = 10, anchor = "w", expand = T, fill = 'both')
tkpack (frame.f, ipady = 10, ipadx = 10, padx = 10, pady = 10, anchor = "w", expand = T, fill = 'both')
tkpack (frame.g, ipady = 10, ipadx = 10, padx = 10, pady = 10, anchor = "w", expand = T, fill = 'both')
tkpack (frame.b, frame.c, side = 'left', ipady = 10, ipadx = 10, padx = 10, pady = 10, expand = T, fill = 'both')
tkpack (frame.d, side = "bottom", pady = 10, padx = 10, expand = T, fill = 'both')

tkpack (spec.frm)
tkbind (base, "<Destroy>", function() tclvalue(cancel)<-2)  

tkraise (base)
tkwait.variable (cancel)
}

