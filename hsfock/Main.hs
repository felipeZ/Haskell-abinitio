
-- The thing...
-- A system shell for the HaskellFock SCF Project 
-- @2012,2013 Angel Alvarez Adhesive tape
-- @2012,2013 Felipe Zapata core SCF machinery 


module Main where

import Data.List (find)
import Data.Maybe ( fromMaybe )
import Control.Concurrent
import Control.Concurrent.Async
import Control.Exception.Base (evaluate)
import Control.Monad.IO.Class
import Control.Monad.Trans.Either
import System.Cmd ( system )
import System.Console.GetOpt
import System.Directory (doesDirectoryExist)
import System.Environment ( getArgs )
import System.FilePath
import System.IO
import Text.Printf

-- Cabal imports
import Data.Version (showVersion)
import Distribution.Version
import Paths_HartreeFock as HsApp

import Options                           -- Here we define the needed fields
import OptsCheck            

import Derivatives                       -- Hartree-Fock Energy derivatives
import HartreeFock                       -- This is the main HartreeFock library, providing functions and types needed
import LinearAlgebra
import Logger                            -- Logger functions 
import Project                           -- We define a generic project skeleton for sample data
import SampleProjects                    -- We include several sample projects for ease of testing


progName = "HaskellFock SCF project , Version:" ++ currVersion ++ " "
    where
        currVersion :: String
        currVersion = showVersion HsApp.version

progAuthors = "@2013 Felipe Zapata, Angel Alvarez, Alessio Valentini"

-- default options
defaultOptions    = Options
 { 
--      optDump        = False
   optModules     = [("scf",doSCF), ("integrals",doIntegrals), ("print",printFiles)]    -- Two payload modules (print is also de default action right now)
 , optMode        = Nothing
 , optVerbose     = False
 , optShowVersion = False
 , optOutput      = Nothing
 , optDataDir     = Nothing
 , optInput       = []
 }

-- currently supported options
acceptedOptions :: [OptsPolicy]
acceptedOptions =
 [ 
   Option ['h','?'] ["help"]    (NoArg  ( check_help           ))                "Show this help message."
 , Option ['V']     ["Version"] (NoArg  ( check_version        ))                "Show version number"
 , Option ['D']     ["datadir"] (ReqArg ( check_data_dir       ) "Dir")          "Directory where files are located"
 , Option ['m']     ["mode"]    (ReqArg ( check_operation_mode ) "Mode")         "Mode of Operation"
 , Option []        ["dump"]    (NoArg  ( check_dump_options   ))                "Force args cmdline dump"
--  , Option ['e']     ["error"]   (NoArg (\ _opts -> return $ Left "forced error on args detected!"))  "Force args checking error"
 , Option ['v']     ["verbose"] (NoArg  ( check_verbosity      ))                "Verbose run on stderr"
--  , Option ['i']     ["input"]   (OptArg (\f opts -> check_input_file f opts) "FILE")             "Input file"
 ]


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
--     let projectdata = project "CO" "STO-3G"
--         charge = 0
--         atom1 = AtomData r1 baseO  8.0
--         atom2 = AtomData r2 baseC  6.0
--         [r1, r2] = atomList projectdata
--         [baseO,baseC] = pBasis projectdata
--     putStrLn "Starting main SCF calculations, please wait...."
--     putStrLn "number of shells"
--     print .  sum . fmap (length . getBasis) $ [atom1,atom2]
--     putStrLn "overlap"
--     overlap <- mtxOverlap [atom1,atom2]
--     print overlap
--     print "Hcore"
--     core <- hcore [atom1,atom2]
--     print core
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

doSCF :: Options -> IO ()
doSCF optsR = do
    logger <- initLogger "water_sto_6-31Gs.out"
    let projectdata =project "water" "6-31G*" 
        charge = 0
        atom1 = AtomData r1 baseO 8.0
        atom2 = AtomData r2 baseH 1.0
        atom3 = AtomData r3 baseH 1.0
        [r1, r2, r3] = atomList projectdata
        [baseH,baseO] = pBasis projectdata
    logMessage logger "Starting main SCF calculations, please wait....\n"
    logMessage logger "number of shells: "
    logMessage logger $ printf "%d\n" $  sum . fmap (length . getBasis) $ [atom1,atom2,atom3]
    hartreeData <- scfHF [atom1,atom2,atom3] charge $ logMessage logger
    logMessage logger "Hartree-Fock has succeeded !!!\n"
    logMessage logger "HF\n"
    logMessage logger $ printf "%.8f\n" $ getEnergy hartreeData
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
--     hartreeData <- scfHF [atom1,atom2,atom3] charge $ logMessage logger
--     logMessage logger "Hartree-Fock has succeeded !!!\n"
--     logMessage logger "HF\n"
--     logMessage logger $ printf "%.8f\n" $ getEnergy hartreeData
--     -- logMessage logger "Calculating the gradient\n"
--     -- gradient <- energyGradient [atom1,atom2,atom3] hartreeData
--     -- logMessage logger $ show gradient
--     logStop logger
-- --     print $ fmap (\x -> checknonZeroDerivative [atom1,atom2,atom2] x [6,6]) [0..8] 

doIntegrals :: Options -> IO ()
doIntegrals optsR = do
    logger <- initLogger "Integrals_water_sto_3g.out"
    let projectdata = project "water" "STO-3G"
        charge = 0
        atom1 = AtomData r1 baseO 8.0
        atom2 = AtomData r2 baseH 1.0
        atom3 = AtomData r3 baseH 1.0
        [r1, r2, r3] = atomList projectdata
        [baseH,baseO] = pBasis projectdata
    logMessage logger "Starting main SCF calculations, please wait....\n"
    logMessage logger "number of shells: "
    logMessage logger $ printf "%d\n" $  sum . fmap (length . getBasis) $ [atom1,atom2,atom3]
    integrals   <- calcIntegrals [atom1,atom2,atom3]
    evaluate integrals
    logMessage logger "DONE"
    logStop logger


-- ============================================= No user serviceable parts beyond this point! ====================================================

main :: IO ()
main = do
    args   <- getArgs
    cores  <- getNumCapabilities
    progHeader cores
    result <- runEitherT $ progOpts args defaultOptions acceptedOptions  -- thread Options over monadic checkers using EitherT over IO
    either somethingIsWrong doSomeStuff result                           --

-- | somethingIsWrong bang in forn tof the user and let him/her hiow to do proper things
somethingIsWrong :: String -> IO ()    
somethingIsWrong msg = do
             putStrLn $ "\nError: " ++ msg ++ "\n"
             putStrLn $ usageInfo header acceptedOptions

-- | do something after unpacking desired lambda from the options record
doSomeStuff :: Options -> IO ()
doSomeStuff optsR@Options { optMode = mode } = do
    case mode of
         Nothing -> printFiles optsR -- | We can also call doNothing here if we dont like printing the files.... 
         Just fun -> fun optsR
         
-- | default dummy action to print filenames after option processing
printFiles :: Options -> IO ()
printFiles opts@Options { optInput = files, optDataDir = datadir } = do
    mapM_ printargs filepaths 
    where
            dir = fromMaybe "./" datadir
            filepaths = zipWith (combine) (cycle [dir]) files
            printargs :: String -> IO ()
            printargs path = putStrLn $ "Processing path: " ++ path ++ "..."

-- | Silly dummy function in case no payload modele was selected
doNothing :: Options -> IO ()
doNothing _ = somethingIsWrong "There is not payload module to execute!\n"             

-- | Keep calm and curry on, we are the good guys....
progHeader :: Int -> IO ()
progHeader c = 
    putStrLn $ progName ++" " ++ progAuthors ++ "\n" ++ show(c) ++ " processor " ++ (core2string c) ++ " detected."
    where
        core2string :: Int -> String
        core2string c = case c > 1 of
                             True -> "cores"
                             False -> "core"

header :: String
header = "Usage: Options [OPTION...] files..."


-- =============================================== Monadic Options checkers =======================================
-- getOpt will partially apply against the supplied argument do we can just come over the options record
-- check user wants verbosity level increased

-- | check_verbosity activates verbose mode of operation
check_verbosity :: Options -> OptsResult
check_verbosity optsR = do
    liftIO $ putStrLn "Verbose mode selected\n"
    return $ optsR { optVerbose = True }          
          
-- | check_version activates versioning info
check_version :: Options -> OptsResult
check_version optsR = return $ optsR { optShowVersion = True }
          
-- | force help message generation, dont panic we are getting into a left cause, we are just not going anywhere after showing the help
check_help :: Options -> OptsResult
check_help _ = left "Command line help, requested"
    

-- | User passed some directory, we make sure this dir exits so file will matched against it
check_data_dir :: String -> Options -> OptsResult
check_data_dir dir optsR = do
    test <- liftIO $ doesDirectoryExist dir
    case test of
         True -> return $ optsR { optDataDir =  Just dir } 
         False -> left $ "Data directory " ++ dir ++ " does not exist" 


-- | check mode of operation. A list of modules is provided in the options record
check_operation_mode :: String -> Options -> OptsResult 
check_operation_mode mode optsR@Options { optModules = modules } = return $ optsR { optMode = selectedModule }
    where
        selectedModule = case (findmodule mode modules) of
                              Just (_,fun) -> Just fun
                              Nothing      -> Nothing
        findmodule :: String -> [(String, (Options-> IO()))] -> Maybe (String,(Options -> IO ()))
        findmodule mode = find ((== mode).fst)  


-- | dump either options or errors 
check_dump_options :: Options -> OptsResult
check_dump_options optsR = do
    liftIO $ putStrLn $ "\n\n\t" ++ show optsR ++ "\n"
    left "Options record dump, requested."
    
    
