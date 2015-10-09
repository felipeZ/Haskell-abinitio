{-# LANGUAGE DeriveDataTypeable, RecordWildCards #-}

-- A system shell for the HaskellFock SCF Project 
-- @2012-2015 Angel Alvarez Adhesive tape
-- @2012-2015 Felipe Zapata core SCF machinery 


module Main where

-- =============================> Standard and third party libraries <===============================
import Control.Monad (when)
import Control.Concurrent (getNumCapabilities)
import System.Console.CmdArgs
import System.Directory ( doesDirectoryExist, doesFileExist )
import System.Environment (getArgs, withArgs)
import System.FilePath
import System.IO
import Text.Printf


-- Cabal imports
import Data.Version (showVersion)
import Distribution.Version
import Paths_HartreeFock as HsApp

-- =================> Internal Modules <======================
import Science.QuantumChemistry.ConcurrencyTools.Logger    -- Logger functions 
import Science.QuantumChemistry.GlobalTypes

-- import Science.QuantumChemistry.HartreeFock.Derivatives    -- Hartree-Fock Energy derivatives
import Science.QuantumChemistry.HartreeFock.HartreeFock    -- This is the main HartreeFock library, providing functions and types needed
-- import Science.QuantumChemistry.HsFock.OptsCheck            

import Science.QuantumChemistry.HsFock.Project             -- We define a generic project skeleton for sample data
import Science.QuantumChemistry.HsFock.SampleProjects      -- We include several sample projects for ease of testing


-- ============================================================
data HSFOCK =  HSFOCK {
                    scf       :: Bool
                   ,basis     :: String
                   ,xyz       :: FilePath
                   ,outFile   :: FilePath
                   -- ,inputFile :: FilePath

                   } deriving (Show, Data, Typeable)

hsFock = HSFOCK
                  {scf   = def &= help "Solves The HF equation using the Self Consistent Field (SCF) Method" 
                  ,basis = def &= help "Uses basis set, for instance: 6-31G*" &= typ
"BASIS"
                  ,xyz   = def &= help "Use the Cartesian coordinates of the file" &= typFile
                  ,outFile = def &= help "Output file" &= typFile                            
                  }
                  &= summary (progName ++ progAuthors) 
                  &= help "The Hartree-Fock method implemented in Haskell"

progName = printf "HaskellAbInitio v%s\n" currVersion
    where
        currVersion :: String
        currVersion = showVersion HsApp.version

progAuthors = "@2015 Felipe Zapata, Angel Alvarez, Alessio Valentini"

-- | Keep calm and curry on
processors :: Int -> IO ()
processors c = 
    printf "%d processor %s detected.\n" (show c) (core2string c)
    where
        core2string :: Int -> String
        core2string c =if  c > 1 then "cores" else "core"

-- | Sanity check of arguments
checkOpts :: HSFOCK -> IO ()
checkOpts HSFOCK{..}= do
  when (not scf)    $ putStrLn "We only Serve HartreeFock SCF (for the moment)"
  when (null basis) $ putStrLn "You must provide a basis set name"
  --todo check that basis set exist
  [xyzBool,outBool]   <- mapM doesFileExist [xyz,outFile]
  when (xyzBool)  $ printf "Coordinates File %s does not exist" xyz 
  when (outBool)  $ printf "file %s is going to be overwritten\n" outFile 

-- ========================> <========================================

main :: IO ()
main  = do
  args <- getArgs
  -- If the user did not specify any arguments, pretend as "--help" was given
  opts <- (if null args then withArgs ["--help"] else id) (cmdArgs hsFock)
  checkOpts opts
  cores  <- getNumCapabilities
  processors cores
  print opts




-- Main routine for SCF calculations
-- doSCF :: Options -> IO ()
-- doSCF optsR = do
--     let projectdata = project "HHe+" "STO-3G"
--         atom1 = AtomData r1 baseHe 2.0
--         atom2 = AtomData r2 baseH  1.0
--         [r1, r2] = atomList projectdata
--         [baseHe,baseH] = pBasis projectdata
--         charge = 1
--     putStrLn "Starting main SCF calculations, please wait...."
--     putStrLn "number of shells"
--     print .  sum . fmap (length . getBasis) $ [atom1,atom2]
--     putStrLn "overlap"
--     overlap <- mtxOverlap [atom1,atom2]
--     print overlap
--     print "Hcore"
--     core <- hcore [atom1,atom2]
--     print core
--     print "Integrals"
--     integrals <- calcIntegrals [atom1,atom2]
--     print integrals
--     result <- scfHF [atom1,atom2] charge
--     print "HF"
--     print $ getEnergy result
--     print "Final Density"
--     print $ getDensity result
--     print "Orbital Energy"
--     print $ getOrbE result
--     print "Final Coeeficients"
--     print $ getCoeff result
--     print "Final Fock Matrix"
--     print $ getFock result


-- doSCF :: Options -> IO ()
-- doSCF optsR = do
--     logger <- initLogger "CO_sto_3g.out"
--     let projectdata = project "CO" "STO-3G"
--         charge = 0
--         atom1 = AtomData r1 baseO  8.0
--         atom2 = AtomData r2 baseC  6.0
--         [r1, r2] = atomList projectdata
--         [baseO,baseC] = pBasis projectdata
--     logMessage logger "Starting main SCF calculations, please wait....\n"
--     logMessage logger "number of shells: "
--     logMessage logger $ printf "%d\n" $  sum . fmap (length . getBasis) $ [atom1,atom2]
--     hartreeData <- scfHF [atom1,atom2] charge $ logMessage logger
--     logMessage logger "Hartree-Fock has succeeded !!!\n"
--     logMessage logger "HF\n"
--     logMessage logger $ printf "%.8f\n" $ getEnergy hartreeData
--     -- logMessage logger "Calculating the gradient\n"
--     -- gradient <- energyGradient [atom1,atom2,atom3] hartreeData
--     -- logMessage logger $ show gradient
--     logMessage logger "The answer is 42!!"
--     logStop logger


doSCF :: IO ()
doSCF = do
    logger <- initLogger "water_6-31G*.out"  
    let projectdata = project "water" "6-31G*"
        charge = 0
        atom1 = AtomData r1 baseO 8.0
        atom2 = AtomData r2 baseH 1.0
        atom3 = AtomData r3 baseH 1.0
        [r1, r2, r3] = atomList projectdata
        [baseH,baseO] = pBasis projectdata
    putStrLn "Starting main SCF calculations, please wait...."
    logMessage logger "number of shells: "
    logMessage logger $ printf "%d\n" $  sum . fmap (length . getBasis) $ [atom1,atom2,atom3]
    hartreeData <- scfHF [atom1,atom2,atom3] charge $ logMessage logger
    putStrLn "DONE"
    logMessage logger "Hartree-Fock has succeeded !!!\n"
    logMessage logger "HF\n"
    logMessage logger $ printf "%.8f\n" $ getEnergy hartreeData
    logMessage logger "The answer is 42!!"
    logStop logger


-- doSCF :: Options -> IO ()
-- doSCF optsR = do
--     logger <- initLogger "water_sto_3g.out"
--     let projectdata = project "water" "STO-3G"
--         charge = 0
--         atom1 = AtomData r1 baseO 8.0
--         atom2 = AtomData r2 baseH 1.0
--         atom3 = AtomData r3 baseH 1.0
--         [r1, r2, r3] = atomList projectdata
--         [baseH,baseO] = pBasis projectdata
--     logMessage logger "Starting main SCF calculations, please wait....\n"
--     logMessage logger "number of shells: "
--     logMessage logger $ printf "%d\n" $  sum . fmap (length . getBasis) $ [atom1,atom2,atom3]
--     mtxS <- mtxOverlap  [atom1,atom2,atom3]
--     logMessage logger $ show mtxS
--     hartreeData <- scfHF [atom1,atom2,atom3] charge $ logMessage logger
--     logMessage logger "Hartree-Fock has succeeded !!!\n"
--     logMessage logger "HF\n"
--     logMessage logger $ printf "%.8f\n" $ getEnergy hartreeData
--     -- logMessage logger "Calculating the gradient\n"
--     -- gradient <- energyGradient [atom1,atom2,atom3] hartreeData
--     -- logMessage logger $ show gradient
--     logMessage logger "The answer is 42!!"
--     logStop logger    


