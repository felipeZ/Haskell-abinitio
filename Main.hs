{-# LANGUAGE RecordWildCards #-}

{-|
Module: Main
Description: Initialization of the simulation
Copyright: @2013 Felipe Zapata, Angel Alvarez
           @2016 Felipe Zapata
The HaskellFock SCF Project
A system shell for the HaskellFock SCF Project 
-}


module Main where

-- ===================> Standard and third party libraries <===================
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

import Science.QuantumChemistry.Integrals.IntegralsEvaluation  (hcore)
import Science.QuantumChemistry.NumericalTools.TableBoys   (generateGridBoys)
-- ============================================================

-- | Options for the mode to execute simulation
hsFock = HSFOCK
                  {scf   = def &= help "Solves The HF equation using the Self Consistent Field (SCF) Method" 
                  ,basis = def &= help "Uses basis set, for instance: 6-31G*" &= typ "BASIS"
                  ,xyz   = def &= help "Use the Cartesian coordinates of the file" &= typFile
                  ,charge= 0   &= help "Molecular charge" &= typ "TotalCharge"
                  ,multi = 1   &= help "Spin state" &= typ "SpinState"
                  ,ndiis = 5   &= help "Number of previuos step to store for DIIS" &= typ "Int"
                  ,outFile = "scf.out" &= help "Output file" &= typFile
                  }

-- | Options for the mode configuring the basis set
basisConfig = BasisConfig
   { basisPath = def &= help "Path to the File that containts the basis set format as plain text" &= typFile
     }

{-| There are currently two modes of execute HSFOCK. the first want does a Self-consistent filed to compute
the energy of a molecule. it takes as argument the name of the basis set to use, the path to the molecular
geometry in xyz and an optional output name.

The second mode of operation transform the basis set written in plain text to binary representation for
fast retrieval. This step should only be run once during the configuration of the library.
-}
hsModes :: Mode (CmdArgs HSFOCK)
hsModes = cmdArgsMode $ modes [hsFock, basisConfig]
    &= verbosityArgs [explicit, name "Verbose", name "V"] []
    &= versionArg [explicit, name "version", name "v", summary "The Hartree-Fock method implemented in Haskell"]
    &= summary (progName ++ progAuthors) 
    &= help "The Hartree-Fock method implemented in Haskell"
    &= helpArg [explicit, name "help", name "h"]
    &= program "HsFock"
 
   
-- ==========================<>=============================

progName = printf "HaskellAbInitio v%s\n" currVersion
    where
        currVersion :: String
        currVersion = showVersion HsApp.version

progAuthors = "@2012-2016 Felipe Zapata, Angel Alvarez, Alessio Valentini"

-- | Check number of processors
processors :: Int -> IO ()
processors c = 
    printf "%s processor %s detected.\n" (show c) (core2string c)
    where
        core2string :: Int -> String
        core2string c =if  c > 1 then "cores" else "core"


-- ========================> <========================================

main :: IO ()
main  = do
  args <- getArgs
  -- If the user did not specify any arguments, pretend as if "--help" was given
  opts <- (if null args then withArgs ["--help"] else id) (cmdArgsRun hsModes)
  checkOpts opts
  cores  <- getNumCapabilities
  processors cores
  executeMode opts
  
-- | Execute one of the execution modes
executeMode :: HSFOCK -> IO ()
executeMode opts@HSFOCK{..}      = doSCF opts
executeMode opts@BasisConfig{..} = doBasisConfig opts

-- | Run the actual SCF using the Hartree-Fock calculation
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
