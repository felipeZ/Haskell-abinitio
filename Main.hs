{-# LANGUAGE DeriveDataTypeable, RecordWildCards #-}

-- A system shell for the HaskellFock SCF Project 
-- @2012-2015 Angel Alvarez Adhesive tape
-- @2012-2015 Felipe Zapata core SCF machinery 


module Main where

-- =============================> Standard and third party libraries <===============================
import Control.Monad (when)
import Control.Concurrent (getNumCapabilities)
import System.Console.CmdArgs
import System.Directory (getHomeDirectory)
import System.Environment (getArgs, withArgs)
import Text.Printf


-- Cabal imports
import Data.Version (showVersion)
import Distribution.Version
import Paths_HartreeFock as HsApp

-- =================> Internal Modules <======================
import Science.QuantumChemistry.BasisSet.SerializeBasis    -- process a plain text basis set to binary
import Science.QuantumChemistry.ConcurrencyTools.Logger    -- Logger functions 
import Science.QuantumChemistry.GlobalTypes

-- import Science.QuantumChemistry.HartreeFock.Derivatives    -- Hartree-Fock Energy derivatives
import Science.QuantumChemistry.HartreeFock.HartreeFock    -- This is the main HartreeFock library, providing functions and types needed
-- import Science.QuantumChemistry.HsFock.OptsCheck            

import Science.QuantumChemistry.HsFock.Initialize          -- Initialize the Atom data type
import Science.QuantumChemistry.HsFock.OptsCheck           -- Sanity options check

import Science.QuantumChemistry.HsFock.Project             -- We define a generic project skeleton for sample data
import Science.QuantumChemistry.HsFock.SampleProjects      -- We include several sample projects for ease of testing


-- ============================================================


hsFock = HSFOCK
                  {scf   = def &= help "Solves The HF equation using the Self Consistent Field (SCF) Method" 
                  ,basis = def &= help "Uses basis set, for instance: 6-31G*" &= typ
"BASIS"
                  ,xyz   = def &= help "Use the Cartesian coordinates of the file" &= typFile
                  ,charge= 0   &= help "Molecular charge" &= typ "TotalCharge"
                  ,multi = 1   &= help "Spin state" &= typ "SpinState"
                  ,outFile = "scf.out" &= help "Output file" &= typFile
                  }
                  &= summary (progName ++ progAuthors) 
                  &= help "The Hartree-Fock method implemented in Haskell"


basisConfig = BasisConfig
   { basisPath = def &= help "Path to the File that containts the basis set format as plain text"
     }
  
  
-- =======================<>=================================================

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


-- ========================> <========================================

main :: IO ()
main  = do
  args <- getArgs
  -- If the user did not specify any arguments, pretend as "--help" was given
  opts <- (if null args then withArgs ["--help"] else id) (cmdArgs hsFock)
  checkOpts opts
  cores  <- getNumCapabilities
  processors cores
  doSomething opts

doSomething :: HSFOCK -> IO ()
doSomething opts@HSFOCK{..}      = doSCF opts
doSomething opts@BasisConfig{..} = doBasisConfig opts

doSCF :: HSFOCK -> IO ()
doSCF hs@HSFOCK{..} = do
  putStrLn "Starting main SCF calculations, please wait...."
  logger      <- initLogger outFile
  atoms       <- initializeAtoms hs (logMessage logger)
  hartreeData <- scfHF atoms charge (logMessage logger)
  logMessage logger "Hartree-Fock has succeeded !!!\n"
  logStop logger


doBasisConfig :: HSFOCK -> IO ()
doBasisConfig opts@BasisConfig{..} = serializeBasisFile basisPath

  -- if (not . null $  basisPath )
  --    then return basisPath
  --    else defaultPathBasis  



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


