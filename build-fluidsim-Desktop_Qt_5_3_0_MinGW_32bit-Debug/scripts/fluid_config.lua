settings = {

  -- constants within sphfluidsimulation
  ratioOfSpecificHeats           = 1.0,
  pressureCoefficient            = 20.0,
  initialDensity                 = 20.0,
  viscosityCoefficient           = 0.03,
  particleMass                   = 1.0,
  maximumVelocity                = 100.0,
  maximumAcceleration            = 100.0,
  motionDampingCoefficient       = 0.0,
  boundaryDampingCoefficient     = 0.2,
  gravityMagnitude               = 9.8,
  isMotionDampingEnabled         = true,
  isBoundaryParticlesEnabled     = false,
  displaySimulationConsoleOutput = true,

  -- config for setting up simulation
  fps                            = 30,
  smoothingRadius                = 0.2,
  numParticles                   = 16000,
  initialDampingConstant         = 8.0,
  finalDampingConstant           = 0.0,
  minColorDensity                = 0.0,
  maxColorDensity                = 100,
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
