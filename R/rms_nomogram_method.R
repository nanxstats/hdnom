# print.nomogram.raw -----------------------------------------------------------

# The following code is derived from the code at
# https://github.com/harrelfe/rms/blob/master/R/nomogram.s
# which is a part of the R package rms released under GPL,
# originally written by Frank E Harrell Jr.

#' @export
print.nomogram.raw <- function(x, dec = 0, ...) {
  obj <- x
  w <- diff(range(obj$lp$x)) / diff(range(obj$lp$x.real))
  cat(
    "Points per unit of linear predictor:", format(w),
    "\nLinear predictor units per point   :", format(1 / w), "\n\n"
  )

  fun <- FALSE
  for (x in names(obj)) {
    k <- x == "total.points" || x == "lp" || x == "abbrev"
    if (k) {
      fun <- TRUE
      next
    }
    y <- obj[[x]]
    if (fun) {
      z <- cbind(round(y[[1]], dec), y$x.real)
      dimnames(z) <- list(rep("", nrow(z)), c("Total Points", x))
    } else {
      z <- cbind(format(y[[1]]), format(round(y$points, dec)))
      dimnames(z) <- list(rep("", length(y$points)), c(x, "Points"))
      # didn't use data.frame since wanted blank row names
    }
    cat("\n")
    print(z, quote = FALSE)
    cat("\n")
  }
  invisible(x)
}

# plot.nomogram.raw ------------------------------------------------------------

# The following code is derived from the code at
# https://github.com/harrelfe/rms/blob/master/R/plot.nomogram.s
# which is a part of the R package rms released under GPL,
# originally written by Frank E Harrell Jr.

#' @importFrom graphics axis lines par segments strwidth text

plot.nomogram.raw <- function(
  x, lplabel = "Linear Predictor",
  fun.side, col.conf = c(1, .3),
  conf.space = c(.08, .2),
  label.every = 1, force.label = FALSE,
  xfrac = .35, cex.axis = .85, cex.var = 1,
  col.grid = NULL,
  varname.label = TRUE, varname.label.sep = "=", ia.space = .7,
  tck = NA, tcl = -0.25, lmgp = .4, naxes,
  points.label = "Points", total.points.label = "Total Points",
  total.sep.page = FALSE, total.fun, cap.labels = FALSE, ...) {
  set <- x

  info <- attr(set, "info")
  fun <- info$fun
  fun.at <- info$fun.at
  nfun <- length(fun)
  funlabel <- info$funlabel
  fun.at <- info$fun.at
  fun.lp.at <- info$fun.lp.at
  R <- info$R
  sc <- info$sc
  maxscale <- info$maxscale
  Intercept <- info$Intercept
  Abbrev <- info$Abbrev
  conf.int <- info$conf.int
  lp <- info$lp
  lp.at <- info$lp.at
  su <- info$space.used
  nint <- info$nint
  discrete <- info$discrete
  minlength <- info$minlength

  col.conf <- rep(col.conf, length = length(conf.int))

  space.used <- su[1] + ia.space * su[2]

  oldpar <- save_oldpar()
  mgp <- oldpar$mgp
  mar <- oldpar$mar
  par(mgp = c(mgp[1], lmgp, mgp[3]), mar = c(mar[1], 1.1, mar[3], mar[4]))
  on.exit(set_par_nro(oldpar))
  tck2 <- tck / 2
  tcl2 <- tcl / 2
  tck3 <- tck / 3
  tcl3 <- tcl / 3

  se <- FALSE
  if (any(conf.int > 0)) {
    se <- TRUE
    zcrit <- qnorm((conf.int + 1) / 2)
    bar <- function(x, y, zcrit, se, col.conf, nlev = 4) {
      y <- rep(seq(y[1], y[2], length = nlev), length.out = length(x))
      for (j in 1:length(x)) {
        xj <- x[j]
        yj <- y[j]
        W <- c(0, zcrit) * se[j]
        for (i in 1:length(zcrit)) {
          segments(xj - W[i + 1], yj, xj - W[i], yj, col = col.conf[i], lwd = 1)
          segments(xj + W[i + 1], yj, xj + W[i], yj, col = col.conf[i], lwd = 1)
        }
      }
    }
  }

  if (!missing(fun.side)) {
    if (!is.list(fun.side)) fun.side <- rep(list(fun.side), nfun)
    if (any(!(unlist(fun.side) %in% c(1, 3)))) {
      stop("fun.side must contain only the numbers 1 and 3")
    }
  }

  num.lines <- 0

  entities <- 0

  # start <- len <- NULL
  # end <- 0

  # Determine how wide the labels can be
  xl <- -xfrac * maxscale

  if (missing(naxes)) {
    naxes <- if (total.sep.page) {
      max(space.used + 1, nfun + lp + 1)
    } else {
      space.used + 1 + nfun + lp + 1
    }
  }

  newpage <- function(
    naxes, xl, maxscale, cex.var, nint, space.used,
    col.grid, cex.axis, tck, tck2, tcl, tcl2,
    label.every, force.label,
    points = TRUE, points.label = "Points", usr) {
    y <- naxes - 1
    plot(
      0, 0,
      xlim = c(xl, maxscale), ylim = c(0, y),
      type = "n", axes = FALSE, xlab = "", ylab = ""
    )
    if (!missing(usr)) par(usr = usr)
    if (!points) return(y + 1)

    ax <- c(0, maxscale)
    text(xl, y, points.label, adj = 0, cex = cex.var)
    x <- pretty(ax, n = nint)
    dif <- x[2] - x[1]
    x2 <- seq((x[1] + x[2]) / 2, max(x), by = dif)
    x2 <- sort(c(x2 - dif / 4, x2, x2 + dif / 4))
    if (length(col.grid)) {
      segments(x, y, x, y - space.used, col = col.grid[1], lwd = 1)
      segments(x2, y, x2, y - space.used, col = col.grid[-1], lwd = 1)
    }
    axisf(
      3,
      at = x, pos = y, cex = cex.axis, tck = tck, tcl = tcl,
      label.every = label.every,
      force.label = force.label, padj = 0
    )
    axisf(
      3,
      at = x2, labels = FALSE, pos = y,
      tck = tck2, tcl = tcl2, cex = cex.axis
    )
    y
  }

  y <- newpage(
    naxes, xl, maxscale, cex.var, nint, space.used, col.grid,
    cex.axis, tck, tck2, tcl, tcl2,
    label.every = label.every,
    force.label = force.label, points.label = points.label
  )

  i <- 0
  ns <- names(set)

  for (S in set[ns %nin% c("lp", "total.points", funlabel)]) {
    i <- i + 1
    setinfo <- attr(S, "info")
    type <- setinfo$type
    y <- y - (if (type == "continuation") ia.space else 1)
    if (y < -.05) {
      y <- newpage(
        naxes, xl, maxscale, cex.var, nint,
        space.used, col.grid,
        cex.axis, tck, tck2, tcl, tcl2,
        label.every = label.every, force.label = force.label,
        points.label = points.label
      ) - (if (type == "continuation") ia.space else 1)
    }
    # word wrap the labels so that they fit into the supplied space.
    text(
      xl, y,
      paste(str_graph_wrap(ns[[i]], abs(xl), cex = cex.var), collapse = "\n"),
      adj = 0, cex = cex.var
    )
    x <- S[[1]]
    nam <- names(S)[1] # stored with fastest first
    if (is.character(x) && nam %in% names(Abbrev)) {
      transvec <- Abbrev[[nam]]$abbrev
      names(transvec) <- Abbrev[[nam]]$full
      x <- transvec[x]
    }

    fx <- if (is.character(x)) {
      x
    } else {
      sedit(formati(x), " ", "")
    } # axis not like bl - was translate()

    xt <- S$points
    # Find flat pieces and combine their labels
    r <- rle(xt)
    if (any(r$length > 1)) {
      is <- 1
      for (j in r$length) {
        ie <- is + j - 1
        if (j > 1) {
          fx[ie] <- if (discrete[nam] || ie < length(xt)) {
            paste(fx[is], "-", fx[ie], sep = "")
          } else {
            paste(fx[is], "+", sep = "")
          }

          fx[is:(ie - 1)] <- ""
          xt[is:(ie - 1)] <- NA
        }
        is <- ie + 1
      }
      fx <- fx[!is.na(xt)]
      xt <- xt[!is.na(xt)]
    }

    # record the side changes
    side <- c(1, 3)
    # subtract 0.6 from the side 1 mgp so that the labels are equaly seperated from the axis
    padj <- c(1, 0)
    new.mgp <- vector(mode = "list", 2)
    new.mgp[[2]] <- c(0, lmgp, 0)
    new.mgp[[1]] <- new.mgp[[2]] - c(0, 0.6, 0)

    # Find direction changes
    ch <- if (length(xt) > 2) {
      c(FALSE, FALSE, diff(diff(xt) > 0) != 0)
    } else {
      rep(FALSE, length(xt))
    }
    if (discrete[nam] && length(xt) > 1) {
      # categorical - alternate adjacent levels
      j <- order(xt)
      lines(range(xt), rep(y, 2)) # make sure all ticks are connected
      for (k in 1:2) {
        is <- j[seq(k, length(j), by = 2)]
        new.labs <- if (cap.labels) capitalize(fx[is]) else fx[is]
        axisf(
          side[k],
          at = xt[is], labels = new.labs, pos = y,
          cex = cex.axis, tck = tck, tcl = tcl,
          force.label = force.label ||
            (minlength == 1 && nam %in% names(Abbrev)),
          disc = TRUE, mgp = new.mgp[[k]], padj = padj[k]
        )
        if (se) {
          bar(
            xt[is],
            if (k == 1) y - conf.space - .32 else y + conf.space + .32,
            zcrit, sc * S$se.fit[is], col.conf
          )
        }
      }
    } else if (!any(ch)) {
      axisf(
        1,
        at = xt, labels = fx, pos = y, cex = cex.axis,
        tck = tck, tcl = tcl, mgp = new.mgp[[1]],
        label.every = label.every, force.label = force.label,
        disc = discrete[nam], padj = padj[1]
      )
      if (se) bar(xt, y + conf.space, zcrit, sc * S$se.fit, col.conf)
    } else {
      lines(range(xt), rep(y, 2)) # make sure all ticks are connected
      j <- (1:length(ch))[ch]
      if (max(j) < length(ch)) j <- c(j, length(ch) + 1)
      flag <- 1
      is <- 1

      for (k in j) {
        ie <- k - 1
        axisf(
          side[flag],
          at = xt[is:ie], labels = fx[is:ie],
          pos = y, cex = cex.axis, tck = tck, tcl = tcl,
          label.every = label.every, force.label = force.label,
          mgp = new.mgp[[flag]],
          disc = discrete[nam], padj = padj[flag]
        )
        if (se) {
          bar(
            xt[is:ie],
            if (side[flag] == 1) {
              y - conf.space - .32
            } else {
              y + conf.space + .32
            },
            zcrit, sc * S$se.fit[is:ie], col.conf
          )
        }
        flag <- if (flag == 2) 1 else 2
        is <- ie + 1
      }
    }
  }

  S <- set$total.points
  x <- S$x
  new.max <- max(x)
  xl.old <- xl
  xl <- -xfrac * new.max
  u <- par()$usr
  if (!missing(total.fun)) total.fun()
  usr <- c(xl * u[1] / xl.old, new.max * u[2] / maxscale, u[3:4])
  par(usr = usr)
  x.double <- seq(x[1], new.max, by = (x[2] - x[1]) / 5)

  y <- y - 1

  if (y < -.05 || total.sep.page) {
    y <- newpage(
      naxes, xl, maxscale, cex.var, nint, space.used, col.grid,
      cex.axis, tck, tck2, tcl, tcl2,
      label.every = label.every, force.label = force.label,
      points = FALSE, usr = usr
    ) - 1
  }

  text(xl, y, total.points.label, adj = 0, cex = cex.var)
  axisf(
    1,
    at = x, pos = y, cex = cex.axis, tck = tck, tcl = tcl,
    label.every = label.every,
    force.label = force.label, mgp = c(0, lmgp - 0.6, 0), padj = 1
  )
  axisf(
    1,
    at = x.double, labels = FALSE, pos = y, tck = tck2,
    tcl = tcl2, cex = cex.axis
  )

  if (lp) {
    S <- set$lp
    x <- S$x
    x2 <- seq(lp.at[1], max(lp.at), by = (lp.at[2] - lp.at[1]) / 2)
    scaled.x2 <- (x2 - Intercept) * sc
    y <- y - 1
    if (y < -.05) {
      y <- newpage(
        naxes, xl, maxscale, cex.var, nint,
        space.used, col.grid,
        cex.axis, tck, tck2, tcl, tcl2,
        label.every = label.every, force.label = force.label,
        points = FALSE, usr = usr
      ) - 1
    }

    text(xl, y, lplabel, adj = 0, cex = cex.var)
    axisf(1,
          at = x, labels = formati(lp.at), pos = y,
          cex = cex.axis, tck = tck, tcl = tcl,
          label.every = label.every, force.label = force.label,
          mgp = c(0, lmgp - 0.6, 0), padj = 1
    )
    axisf(1,
          at = scaled.x2, labels = FALSE, tck = tck2, tcl = tcl2,
          pos = y, cex = cex.axis
    )
    conf <- S$conf
    if (length(conf)) {
      bar(conf$x,
          y + c(conf.space[1], conf.space[1] + conf$w * diff(conf.space)),
          zcrit, conf$se, col.conf,
          nlev = conf$nlev
      )
    }
  }

  i <- 0
  if (nfun > 0) {
    for (S in set[funlabel]) {
      i <- i + 1
      y <- y - 1
      scaled <- S$x
      fat <- S$fat
      s <- S$which # ???
      if (y < -.05) {
        y <- newpage(
          naxes, xl, maxscale, cex.var, nint, space.used,
          col.grid, cex.axis, tck, tck2, tcl, tcl2,
          label.every = label.every, force.label = force.label,
          points = FALSE, usr = usr
        ) - 1
      }
      text(xl, y, funlabel[i], adj = 0, cex = cex.var)
      sides <- if (missing(fun.side)) {
        rep(1, length(fat))
      } else {
        (fun.side[[i]])[s]
      }
      if (length(sides) != length(fat)) {
        stop("fun.side vector not same length as fun.at or fun.lp.at")
      }
      for (jj in 1:length(fat))
        axis(
          sides[jj],
          at = scaled[jj], labels = fat[jj],
          pos = y, cex.axis = cex.axis, tck = tck, tcl = tcl,
          mgp = if (sides[jj] == 1) c(0, lmgp - 0.6, 0) else c(0, lmgp, 0),
          padj = if (sides[jj] == 1) 1 else 0
        )
      lines(range(scaled), rep(y, 2)) # make sure all ticks are connected
    }
  }
  invisible(x)
}

# capitalize the first letter of a string
capitalize <- function(string) {
  capped <- grep("^[A-Z]", string, invert = TRUE)
  substr(string[capped], 1, 1) <- toupper(substr(string[capped], 1, 1))
  string
}

# Version of axis allowing tick mark labels to be forced, or to
# label every 'label.every' tick marks

axisf <- function(
  side, at, labels = TRUE, pos, cex, tck, tcl,
  label.every = 1, force.label = FALSE, disc = FALSE, ...) {
  ax <- function(..., cex) axis(..., cex.axis = cex)

  ax(side, at, labels = FALSE, pos = pos, cex = cex, tck = tck, tcl = tcl, ...)

  if (is.logical(labels) && !labels) return(invisible())

  if (label.every > 1 && !disc) {
    sq <- seq(along = at, by = label.every)
    at[-sq] <- NA
  }
  if (is.logical(labels)) labels <- format(at, trim = TRUE)

  if (force.label) {
    for (i in 1:length(labels))
      if (!is.na(at[i])) {
        ax(side, at[i], labels[i], pos = pos, cex = cex, tcl = 0, ...)
      }
  }
  else {
    ax(side, at[!is.na(at)], labels[!is.na(at)],
       pos = pos, cex = cex, tcl = 0, ...
    )
  }

  invisible()
}

save_oldpar <- function() {
  # saves existing state of par() and makes changes suitable
  # for restoring at the end of a high-level graphics functions
  oldpar <- par()
  oldpar$fin <- NULL
  oldpar$new <- FALSE
  invisible(oldpar)
}

str_graph_wrap <- function(
  x, width = 0.9 * getOption("width"),
  indent = 0, exdent = 0,
  prefix = "", simplify = TRUE, units = "user", cex = NULL) {
  if (!is.character(x)) x <- as.character(x)

  spc.len <- strwidth(" ", units = units, cex = cex)
  prefix.len <- strwidth(prefix, units = units, cex = cex)
  indentString <- paste(rep.int(" ", indent), collapse = "")
  indent <- indent * spc.len
  exdentString <- paste(rep.int(" ", exdent), collapse = "")
  exdent <- exdent * spc.len

  y <- list()
  z <- lapply(strsplit(x, "\n[ \t\n]*\n"), strsplit, "[ \t\n]")
  for (i in seq_along(z)) {
    yi <- character(0)
    for (j in seq_along(z[[i]])) {
      words <- z[[i]][[j]]
      nc <- strwidth(words, units = units, cex = cex)
      if (any(is.na(nc))) {
        nc0 <- strwidth(words, units = units, cex = cex)
        nc[is.na(nc)] <- nc0[is.na(nc)]
      }
      if (any(nc == 0)) {
        zLenInd <- which(nc == 0)
        zLenInd <- zLenInd[!(zLenInd %in% (grep("\\.$", words) + 1))]
        if (length(zLenInd) > 0) {
          words <- words[-zLenInd]
          nc <- nc[-zLenInd]
        }
      }
      if (length(words) == 0) {
        yi <- c(yi, "", prefix)
        next
      }
      currentIndex <- 0
      lowerBlockIndex <- 1
      upperBlockIndex <- integer(0)
      lens <- cumsum(nc + spc.len)
      first <- TRUE
      maxLength <- width - prefix.len - indent
      while (length(lens) > 0) {
        k <- max(sum(lens <= maxLength), 1)
        if (first) {
          first <- FALSE
          maxLength <- maxLength + indent - exdent
        }
        currentIndex <- currentIndex + k
        if (nc[currentIndex] == 0) {
          upperBlockIndex <- c(upperBlockIndex, currentIndex - 1)
        } else {
          upperBlockIndex <- c(upperBlockIndex, currentIndex)
        }
        if (length(lens) > k) {
          if (nc[currentIndex + 1] == 0) {
            currentIndex <- currentIndex + 1
            k <- k + 1
          }
          lowerBlockIndex <- c(lowerBlockIndex, currentIndex + 1)
        }
        if (length(lens) > k) {
          lens <- lens[-(1:k)] - lens[k]
        } else {
          lens <- NULL
        }
      }
      nBlocks <- length(upperBlockIndex)
      s <- paste(prefix, c(indentString, rep.int(
        exdentString,
        nBlocks - 1
      )), sep = "")
      for (k in (1:nBlocks))
        s[k] <- paste(
          s[k], paste(words[lowerBlockIndex[k]:upperBlockIndex[k]], collapse = " "),
          sep = ""
        )
      yi <- c(yi, s, prefix)
    }
    y <- if (length(yi)) {
      c(y, list(yi[-length(yi)]))
    } else {
      c(y, "")
    }
  }
  if (simplify) {
    y <- unlist(y)
  }
  y
}

set_par_nro <- function(pars) {
  # Sets non-read-only par parameters from the input list
  i <- names(pars) %nin%
    c("cin", "cra", "csi", "cxy", "din", "xlog", "ylog", "gamma", "page")
  invisible(par(pars[i]))
}
