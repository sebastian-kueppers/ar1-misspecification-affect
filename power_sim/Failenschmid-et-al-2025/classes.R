#' ----------------------------------------------------------------------------#
#' Title: Classes and Generic Methods                                          #
#' Author: Jan Ian Failenschmid                                                #
#' Created Date: 25-03-2024                                                    #
#' -----                                                                       #
#' Last Modified: 13-02-2025                                                   #
#' Modified By: Jan Ian Failenschmid                                           #
#' -----                                                                       #
#' Copyright (c) 2024 by Jan Ian Failenschmid                                  #
#' E-mail: J.I.Failenschmid@tilburguniveristy.edu                              #
#' -----                                                                       #
#' License: GNU General Public License v3.0 or later                           #
#' License URL: https://www.gnu.org/licenses/gpl-3.0-standalone.html           #
#' ----------------------------------------------------------------------------#



### Define classes, generics and methods ---------------------------------------

#' Generative models
#' The gen_model class formalises and structures generative dynamic time-series
#' models so that they can easily be simulated using sim_tsm. In this
#' simulation it is used to specify the different data generating processes
#' and to simulated from them.

setClass(
  "gen_model",
  slots = list(
    # Model name to be copied into method
    model_name = "character",
    # Model type is either DE: "differential equation"
    # or SSM: "state-space model"
    model_type = "character",
    # End time of the simulation
    time = "numeric",
    # Simulation steps of the Euler-Maruyama Method
    delta = "numeric",
    # Stepsize of the simulated time-series
    stepsize = "numeric",
    # List of parameters of the generative model, these parameters can be used
    # in the following formulas
    pars = "list",
    # List of formulas specifying the initial conditions of the times-series
    start = "list",
    # List of formulas for the deterministic part of the state equation
    state_eq = "list",
    # List of formulas for the stochastic part of the state equation
    dynamic_error = "list",
    # List of formulas for the measurement model
    meas_eq = "list",
    # Function specifying how the boundry conditions should be handled
    boundary = "function"
  )
)

setMethod("initialize", "gen_model", function(
    #' Initializer for objects of the gen_model class
    .Object,
    model_name = character(),
    model_type = character(),
    time = 1,
    delta = numeric(),
    stepsize = 1,
    pars = list(),
    start = list(),
    state_eq = list(),
    dynamic_error = list(),
    meas_eq = list(),
    boundary = function(x) ifelse(is.finite(x), x, NA)) {
  # Name model
  if (length(model_name) == 0) {
    model_name <- "my_model"
    warning("Model name missing model will be named my_model")
  }

  if (model_type == "SSM") {
    delta <- NA_real_
  }

  if (model_type == "DE" && length(stepsize) == 0) {
    stepsize <- delta
  }

  # If dynamic error list is empty, set dynamic error to zero for all variables
  if (length(dynamic_error) == 0 && length(state_eq) > 0) {
    dynamic_error <- lapply(
      lapply(state_eq, "[[", 2),
      function(x) {
        as.formula(paste0(x, " ~ 0"))
      }
    )
  }

  # Match formula order between start, state_eq, and dynamic_error
  start <- start[match(
    sapply(state_eq, "[[", 2),
    sapply(start, "[[", 2)
  )]

  dynamic_error <- dynamic_error[match(
    sapply(state_eq, "[[", 2),
    sapply(dynamic_error, "[[", 2)
  )]

  # Assign slots
  .Object@model_name <- model_name
  .Object@model_type <- model_type
  .Object@time <- time
  .Object@delta <- delta
  .Object@stepsize <- stepsize
  .Object@pars <- pars
  .Object@start <- start
  .Object@state_eq <- state_eq
  .Object@dynamic_error <- dynamic_error
  .Object@meas_eq <- meas_eq
  .Object@boundary <- boundary

  # Validate object
  validObject(.Object)

  # Return
  .Object
})

# Methods
setValidity("gen_model", function(object) {
  #' Validator for objects of the gen_model class

  errors <- character()
  # If model type is not valid, return error
  if (!slot(object, "model_type") %in% c("SSM", "DE")) {
    msg <- c("Model type needs to be either SSM or DE.")
    errors <- c(msg, errors)
  }

  # If model time is not a positive real number
  if (is.na(slot(object, "time")) || slot(object, "time") < 0) {
    msg <- c("Time needs to be a positive real number.")
    errors <- c(msg, errors)
  }

  # If delta is negative or a non-integer devisor of the stepsize, return error
  if (slot(object, "model_type") == "DE") {
    if (is.na(slot(object, "delta")) || slot(object, "delta") <= 0) {
      msg <- c("Delta needs to be a positive real number.")
      errors <- c(msg, errors)
    }
    if (round(
      slot(object, "stepsize") / slot(object, "delta"),
      digits = 10
    ) %% 1 != 0) {
      msg <- paste0(
        "Delta needs to be an integer devisor of stepsize: ",
        slot(object, "stepsize"), ". This can be achieved by setting delta to ",
        slot(object, "stepsize"), "/p."
      )
      errors <- c(msg, errors)
    }
  }

  # If stepsize is non-positive or larger than time return an error
  if (slot(object, "stepsize") <= 0) {
    msg <- c("Stepsize needs to be a positive number smaller than Time.")
    errors <- c(msg, errors)
  }

  # If any dynamic error is negative, return error
  if (!all(mapply(
    function(x, y) {
      start <- eval(y[[3]], slot(object, "pars"))
      names(start) <- deparse(y[[2]])
      eval(x[[3]], list2env(c(start, slot(object, "pars"))))
    },
    x = slot(object, "dynamic_error"), y = slot(object, "start")
  ) >= 0)) {
    msg <- c("All dyamic errors need to be real non-negative values")
    errors <- c(msg, errors)
  }

  if (length(errors) == 0) TRUE else errors
})

setMethod("show", "gen_model", function(object) {
  #' Print method for objects of the gen_model class.

  out <- character()
  if (slot(object, "model_type") == "SSM") {
    out <- c(out, paste0(
      "Generative State Space Model: ",
      slot(object, "model_name"), ".\n"
    ))
    out <- c(out, paste0(
      "The model will be simulated from 0 to ",
      slot(object, "time"), " in steps of ", slot(object, "stepsize"), ".\n"
    ))
  } else {
    out <- c(out, paste0(
      "Generative Differential Equation Model: ",
      slot(object, "model_name"), ".\n"
    ))
    out <- c(out, paste0(
      "The model will be simulated from 0 to ",
      slot(object, "time"), " in steps of ", slot(object, "stepsize"), "."
    ))
    out <- c(out, paste0(
      "The Euler-Maruyama method will use discretization steps of ",
      slot(object, "delta"), "."
    ))
  }

  out <- c(out, paste0(
    "\nState starting values:"
  ))
  out <- c(out, sapply(slot(object, "start"), deparse))

  out <- c(out, paste0(
    "\nDeterministic dynamic state equations:"
  ))
  out <- c(out, sapply(slot(object, "state_eq"), deparse))

  out <- c(out, paste0(
    "\nDynamic error equations:"
  ))
  out <- c(out, sapply(slot(object, "dynamic_error"), deparse))

  out <- c(out, paste0(
    "\nMeasurement Equations:"
  ))
  out <- c(out, sapply(slot(object, "meas_eq"), deparse))

  if (length(slot(object, "pars")) > 0) {
    out <- c(out, paste0(
      "\nModel Parameters:"
    ))

    out <- c(out, paste(
      paste(
        names(slot(object, "pars")),
        slot(object, "pars"),
        sep = " = "
      ),
      sep = " "
    ))
  }

  out_msg <- paste0(out, sep = "\n")
  cat(out_msg)
})

## Analysis methods
#' The analysis methods are classes defined for each analysis method used in
#' the simulation. This is mainly done as a convenience to be able to define
#' the generic "fit" and "calculate_perfomance_measures" functions and assign
#' seperate behaviours via methods to each analysis method. For this, a parent
#' "method" class is created that should never be initialized. Instead it holds
#' the structure for the simulation results and the generic function that can
#' be inherited by each of the specific analysis methods. Each analysis method
#' is then a seperate class that inherits from the method class and has their
#' respective methods specified in a seperate file.

setClass(
  "method",
  slots = list(
    # Method name
    method_name = "character",
    # Name of the respective generative model
    gen_model = "character",
    # Simulation conditions
    conditions = "list",
    # Mean function estimate
    estimate = "numeric",
    # Confidence interval estimate
    ci = "list",
    # Mean squared error estimate
    mse = "numeric",
    # Generalized Cross validation criterion
    gcv = "numeric",
    # Confidence interval coverage
    ci_coverage = "numeric",
    # Wiggliness parameter
    wiggliness = "numeric",
    # Convergence
    converged = "logical"
  ),
  prototype = list(
    method_name = NA_character_,
    gen_model = NA_character_,
    conditions = list(),
    estimate = NA_real_,
    ci = list(lb = NA, ub = NA),
    mse = NA_real_,
    gcv = NA_real_,
    ci_coverage = NA_real_,
    wiggliness = NA_real_,
    converged = FALSE
  )
)

#' Generic function for fitting each of the analyis methods defined for method
#' to be inherited.
setGeneric("fit", function(method, data) standardGeneric("fit"),
  signature = "method"
)

#' Generic function for obtaining performance measures from a fitted analysis
#' method to be inherited.

setMethod("plot", "method", function(
    x, sim = NULL, observation = "y_obs",
    state = "y", ...) {
  #' Plot method defined for method to be inherited by each analysis method.
  #' Can be used as a wrapper to pass the plot call to the fit slot of the
  #' method object or to obtain a custom plot by providing the output
  #' generated by simulate to the sim argument.

  # Exctract data set from sim corresponding to the method x
  data <- sim$dat[which(sapply(sim$method, function(method_list) {
    any(sapply(method_list, identical, y = x))
  }))][[1]]

  plot(x = data$time, data[[observation]], ...)
  lines(
    x = data$time, y = data[[state]],
    lwd = 2.5
  )
  if (slot(x, "converged")) {
    lines(x = data$time, y = slot(x, "estimate"), col = "red", lwd = 2.5)
    lines(
      x = data$time, y = slot(x, "ci")$lb, col = "red", lty = 2,
      lwd = 2.5
    )
    lines(
      x = data$time, y = slot(x, "ci")$ub, col = "red", lty = 2,
      lwd = 2.5
    )
  }
})

setMethod("show", "method", function(object) {
  #' Print method defined for method to be inherited by each analysis method.

  out <- character()
  # Print method name
  out <- c(
    out,
    "Method object for", slot(object, "method_name")
  )
  # Print gen model
  if (is.na(slot(object, "gen_model"))) {
    out <- c(
      out,
      "\n\nNo generative model has been assigned to the method."
    )
  } else {
    out <- c(
      out,
      "\n\nThe generative model for the method is",
      slot(object, "gen_model")
    )
  }
  # Print conditions
  if (length(slot(object, "conditions")) == 0) {
    out <- c(
      out,
      "\n\nNo conditions have been assigned to the method."
    )
  } else {
    out <- c(
      out,
      "\n\nThe simulation conditions for the method are:\n",
      paste(
        paste(
          names(slot(object, "conditions")),
          slot(object, "conditions"),
          sep = " = "
        ),
        sep = " "
      )
    )
  }
  out <- c(
    out,
    "\nMSE:", round(slot(object, "mse"), 3),
    "\nGCV:", round(slot(object, "gcv"), 3),
    "\nCI coverage:", round(slot(object, "ci_coverage"), 3),
    "\n Wiggliness:", round(slot(object, "wiggliness"), 3),
    "\n\nThe state inference is:\n"
  )
  cat(out)
  print(data.frame(
    lb = slot(object, "ci")[["lb"]],
    state = slot(object, "estimate"),
    ub = slot(object, "ci")[["ub"]]
  ))
})
