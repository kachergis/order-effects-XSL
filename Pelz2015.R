# attention switching model from PelzPiantadosiKidd ICDL 2015

# parameters:
# (A) the learning curve that the learner follows as he attends to the stimulus, 
# (B) memory decay rate, as a model of the learner’s short term memory limitations, 
# (C) the cost of switching attention between objects (0.5), and 
# (D) the amount of information that the learner knows about the objects at the start of test. 

# Attention allows accumulation of data about an object, but there is a diminishing return 
# as the learner comes to know everything that can be learned about each object. 
# Importantly, these objects are best thought of as dynamic stimuli that will change when unattended, 
# meaning that information will be lost about them when attention is oriented elsewhere.

# At each timestep, the model computes the amount of information that it expects to gain from each object, 
# and chooses whether to continue to attend to the current object or to switch to another based on the 
# expected information gains

# learning curve: Gompertz
info_gain <- function(t, a=100, c=1) {
  # returns amount of information known about each object at time t, 
  # a is the maximum amount of information that can be learned from each object (fixed at 100 bits), 
  # b characterizes the amount of initial information, 
  # c is the learning rate (fixed at 1.0 for simplicity)
  return( a*exp(-b*exp(-c*t)) )
}

# power law decay of information for unattended objects
decay <- function(t, beta=.1) {
  # t = amount of time since object was last attended
  # beta = information decay parameter
  return( (t+1)^(-beta) )
}

# softmax to determine P(choosing to attend obj i)
softmax <- function(obj_info, chi=30) {
  # chi=0 - random choice, chi=Inf - perfect maximization; 
  # chi=30 corresponds to a 95% chance of choosing an option that giving .1 bits more than its alternative
  exp(chi*obj_info) / sum(exp(chi*obj_info))
}

run_sim <- function(nobj, timesteps=100) {
  obj_info = rep(0, nobj) # object information (could be random starting amount)
  info_traj = matrix(nrow=timesteps, ncol=nobj)
  for(t in 1:timesteps) {
    info_traj[i,]
  }
}

run_sim(10)
