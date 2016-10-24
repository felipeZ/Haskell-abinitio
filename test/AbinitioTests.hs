{-|
Module: Main
Description: Initialization of the simulation
Copyright: @2013 Felipe Zapata, Angel Alvarez
           @2016 Felipe Zapata
Testing The haskell-fock project
-}


import Science.QuantumChemistry.GlobalTypes (AtomData(..))
import Science.QuantumChemistry.HsFock.SampleProjects (project)
import Science.QuantumChemistry.HsFock.Project

main :: IO ()
main = do 
      let projectdata = project "water" "STO-3G"
          charge = 0
          atom1 = AtomData r1 baseO 8.0
          atom2 = AtomData r2 baseH 1.0
          atom3 = AtomData r3 baseH 1.0
          [r1, r2, r3] = atomList projectdata
          [baseH,baseO] = pBasis projectdata
      print baseH

-- doSCF :: HSFOCK -> IO ()
-- doSCF _ = do
--     logger <- initLogger "water_sto_3g.out"
--     let projectdata = project "water" "STO-3G"
--     -- logger <- initLogger "water_6_31G*.out"
--     -- let projectdata = project "water" "6-31G*"
--         charge = 0
--         atom1 = AtomData r1 baseO 8.0
--         atom2 = AtomData r2 baseH 1.0
--         atom3 = AtomData r3 baseH 1.0
--         [r1, r2, r3] = atomList projectdata
--         [baseH,baseO] = pBasis projectdata
--         -- gridBoys   = generateGridBoys 0.1 -- ^gridBoys dx, where mMax is the maximum order of the boys function and dx the grid delta
--     logMessage logger "Starting main SCF calculations, please wait....\n"
--     logMessage logger "Core Matrix: "
--     -- core <- hcore gridBoys [atom1,atom2,atom3]
--     -- mtxS <- mtxOverlap  [atom1,atom2,atom3]
--     -- logMessage logger $ show core
--     hartreeData <- scfHF [atom1,atom2,atom3] charge $ logMessage logger
--     logMessage logger "Hartree-Fock has succeeded !!!\n"
--     logMessage logger "HF\n"
--     logMessage logger $ printf "%.8f\n" $ getEnergy hartreeData
-- --     logMessage logger "Calculating the gradient\n"
-- --     -- gradient <- energyGradient [atom1,atom2,atom3] hartreeData
-- --     -- logMessage logger $ show gradient
-- --     logMessage logger "The answer is 42!!"
-- --     logStop logger    
