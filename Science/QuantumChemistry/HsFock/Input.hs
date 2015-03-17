
module Science.QuantumChemistry.HsFock.Input (
              JobInfo(..)
             ,getInfo
             ) where


-- ######################  HERE IS PARSED THE INPUT DATA #################################

-- The InputParse Module is intended for parsing the input file
-- which contains the keywords and geometry of the molecule to simulate

import Text.ParserCombinators.Parsec
import Data.Char (toLower,toUpper)
import Control.Arrow ((***),second)
import Control.Monad(liftM)
import Data.Maybe (fromMaybe)
import qualified Data.Map as M

-- ------------------------> <---------------------------------
import Science.QuantumChemistry.GlobalTypes

-- =====================> TYPES <============================
    
data JobInfo = JobInfo {
    jobTheory :: String
   ,jobBasis  :: [Basis]
   ,jobCharge :: Int
   ,jobElecs  :: Int
   ,jobCoord  :: [(String,[Double])]
    } deriving Show


data InputError = UnkownAtomicElement 

    
-- =========================> Body <===================

getInfo :: String -> Either ParseError JobInfo
getInfo input = do
      dat <- rawParser input
      let fun = (toks *** readerCoord) . breaker
          ([(theory,t),(basis,b),(charge,ch')],geom) = fun dat
          ch = readInteger ch'
          atomicLabels = fst `fmap` geom
          ztotal = sum . map (\k -> fromMaybe msg $ M.lookup (f2U k) atom2charge) $ atomicLabels
          es =  ztotal - ch
      return $ JobInfo t (readBasis b) ch es geom
      
  where msg = error "Sorry boy, but I don't know some of your input atoms"
        f2U (c:cs) = toUpper c : cs
          
readerCoord :: [String] -> [(String,[Double])]
readerCoord = map ((\(x:xs) -> ( x, fun xs)) . words) . tail
  where fun w = map readDouble w

breaker :: [[String]] -> ([String],[String])
breaker = break ("coordinates"==) . filter (not . null) . (map (map toLower)) . concat

toks :: [String] ->[(String,String)]
toks = map (second tail .  break ('=' ==))

readBasis :: String -> [Basis]
readBasis = const [[CGF sto3gHe S]]
  where sto3gHe = [ (0.182985,0.480844),(0.587136,1.776691),(0.607075,9.753934)]

-- ============> RAW PARSER <===================
rawParser :: String -> Either ParseError [[String]]
rawParser input = parse body "Error in the input" input

body :: GenParser Char st [[String]]
body = sepBy line (char '\n')

line :: GenParser Char st [String]
line = sepBy (many (noneOf "\n")) (char '=')

-- ========================> Auxiliar <==============

readDouble  x = read x :: Double
readInteger x = read x :: Int


-- =======================> 
  
