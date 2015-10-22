

module Science.QuantumChemistry.BasisSet.SerializeBasis where

-- ==================> Standard and third party libraries <=================
import Control.Concurrent.Async 
import qualified Data.ByteString as B
import Data.Serialize
import System.Directory (getDirectoryContents)

-- =================> Internal Modules <======================
import Science.QuantumChemistry.ParsecTools.ParserBasis (parseBasisFile)

-- =======================> <==============================

serializeBasisFile :: FilePath -> IO ()
serializeBasisFile = undefined
 -- do 
 --   names <- getDirectoryContents topdir
 --   let properNames = filter (`notElem` [".", ".."]) names
 --   basis <- mapConcurrently (encode . parseBasisFile) properNames 
 --   writeBasisBin basis
