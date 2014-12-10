-- note: color velocity, fading stuck boundary particles, z sorting

settings = {

  -- constants within sphfluidsimulation
  ratioOfSpecificHeats             = 1.0,
  pressureCoefficient              = 20.0,
  initialDensity                   = 20.0,
  viscosityCoefficient             = 0.03,
  particleMass                     = 1.0,
  maximumVelocity                  = 75.0,
  maximumAcceleration              = 75.0,
  motionDampingCoefficient         = 0.0,
  boundaryDampingCoefficient       = 0.6,
  gravityMagnitude                 = 9.8,
  isMotionDampingEnabled           = true,
  isBoundaryParticlesEnabled       = false,
  displaySimulationConsoleOutput   = true,
  
  -- graphics
  minColorDensity                = 0.0,
  maxColorDensity                = 100,
  maxColorVelocity               = 100.0,
  maxColorAcceleration           = 100.0,
  colorArrivalRadius             = 0.1,
  stuckToBoundaryRadius          = 0.01,
  stuckToBoundaryAlphaVelocity   = 3.0,
  isHiddenBoundaryParticlesEnabled = true,

  -- config for setting up simulation
  logfile                        = "file1.txt",
  fps                            = 30,
  smoothingRadius                = 0.2,
  numParticles                   = 1000,
  initialDampingConstant         = 2.0,
  finalDampingConstant           = 0.0,
  isRenderingEnabled             = true,
  isSimulationPaused             = false,
  
  initialBounds                  = {minx = 0.0,
                                    maxx = 3.0,
                                    miny = 0.0,
                                    maxy = 8.0,
                                    minz = 0.0,
                                    maxz = 2.0},
                                    
  finalBounds                    = {minx = 0.0,
                                    maxx = 8.0,
                                    miny = 0.0,
                                    maxy = 8.0,
                                    minz = 0.0,
                                    maxz = 2.0}
    
}

settings.minColorDensity = settings.initialDensity
