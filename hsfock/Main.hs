
-- The thing...
-- A system shell for the HaskellFock SCF Project 
-- @2012 Angel Alvarez System Shell and Data constructor
-- @2012 Felipe Zapata SCF machinery 

module Main where

-- Main libraries get imported here
import Control.Monad
-- import Control.Monad(foldM,liftM,ap)
-- import Data.Maybe ( fromMaybe )
import Data.Either
import Data.Version (showVersion)
import System.Environment
import System.Console.GetOpt

-- Cabal imports
import Distribution.Version
import Paths_HartreeFock as HsFock

-- Our module imports
import OptsCheck
import HartreeFock

-- Keep calm and curry on, we are the good guys....
progHeader :: IO ()
progHeader = 
    putStrLn $ "HaskellFock SCF project, Version:" ++ currVersion ++ " @2012 Felipe Zapata, Angel Alvarez" ++ "\n"
    where
        currVersion :: String
        currVersion = showVersion HsFock.version

header :: String
header = "Usage: Options [OPTION...] files..."

-- default options
defaultOptions    = Options
 { optDump        = False
 , optVerbose     = False
 , optShowVersion = False
 , optOutput      = Nothing
 , optInput       = Nothing
 , optDataDir     = Nothing
 }

-- curently supported options
acceptedOptions :: [OptDescr (Either String Options -> IO (Either String Options))]
acceptedOptions =
 [ 
   Option ['e']     ["error"]   (NoArg (\ _opts -> return $ Left "forced error on args detected!"))  "Force args checking error"
 , Option ['h','?'] ["help"]    (NoArg (\ opts -> fun_check_help opts))                              "Show this help message."
 , Option ['d']     ["dump"]    (NoArg (\ opts -> fun_check_dump_options opts))                      "Force args cmdline dump"
 , Option ['D']     ["datadir"] (ReqArg (\d opts -> fun_check_data_dir d opts) "DIR")                "Data dir"
 , Option ['v']     ["verbose"] (NoArg (\ opts -> fun_check_verbosity opts))                         "Verbose run on stderr"
 , Option ['V']     ["Version"] (NoArg (\ opts -> fun_check_version opts))                           "Show version number"
{- , Option ['i']     ["input"]   (OptArg (\f opts -> fun_check_input_file f opts) "FILE")             "Input file"-}
 ]

main :: IO ()
main = do
    progHeader
    args <- getArgs
    result <- progOpts args defaultOptions acceptedOptions header
    case result of
         Left msg -> do
             putStrLn $ msg ++ "\n"
             putStrLn $ usageInfo header acceptedOptions
         Right options -> do
             let sto3g = [(0.444635,0.168856),(0.535328,0.623913),(0.154329,3.42525)]
                 sto3gHe = [ (0.182985,0.480844),(0.587136,1.776691),(0.607075,9.753934)]
                 r1 = [0.0,0.0,1.4632]
                 r2 = [0.0,0.0,0.0]
                 base1 = [CGF sto3gHe S]
                 base2 = normaCoeff `fmap` [CGF sto3g S]
                 listBasis = [base1, base2]
                 listCoord = [r1,r2]
                 zlist = [2.0,1.0]
             result <- scfHF listCoord listBasis zlist 2
             print "HF"
             print $ getEnergy result
             print "Final Density"
             print $ getDensity result
             print "Orbital Energy"
             print $ getOrbE result
             print "Final Coeeficients"
             print $ getCoeff result
             print "Final Fock Matrix"
             print $ getFock result