// ***************************************************************************
// BamIndex.h (c) 2009 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 19 November 2010 (DB)
// ---------------------------------------------------------------------------
// Provides basic BAM index interface
// ***************************************************************************

#ifndef BAM_INDEX_H
#define BAM_INDEX_H

#include <api_global.h>
#include <BamAux.h>
#include <iostream>
#include <string>
#include <vector>

namespace BamTools {

class BamReader;
class BgzfData;

namespace Internal {
  class BamStandardIndex;
  class BamToolsIndex;
} // namespace Internal

// --------------------------------------------------  
// BamIndex base class
class API_EXPORT BamIndex {

    // specify index-caching behavior
    //
    // @FullIndexCaching - store entire index file contents in memory
    // @LimitedIndexCaching - store only index data for current reference
    //   being processed
    // @NoIndexCaching - do not store any index data.  Load as needed to 
    //   calculate jump offset
    public: enum BamIndexCacheMode { FullIndexCaching = 0
                                   , LimitedIndexCaching
                                   , NoIndexCaching
                                   };
  
    // ctor & dtor
    public:
        BamIndex(BamTools::BgzfData* bgzf, BamTools::BamReader* reader);
        virtual ~BamIndex(void);
        
    // index interface
    public:
        // creates index data (in-memory) from current reader data
        virtual bool Build(void) =0;
        // returns supported file extension
        virtual const std::string Extension(void) const =0;
        // returns whether reference has alignments or no
        virtual bool HasAlignments(const int& referenceID) const =0;
        // attempts to use index to jump to region; returns success/fail
        // a "successful" jump indicates no error, but not whether this region has data
        //   * thus, the method sets a flag to indicate whether there are alignments 
        //     available after the jump position
        virtual bool Jump(const BamTools::BamRegion& region, bool* hasAlignmentsInRegion) =0;
        // loads existing data from file into memory
        virtual bool Load(const std::string& filename);
        // change the index caching behavior
        virtual void SetCacheMode(const BamIndexCacheMode mode);
        // writes in-memory index data out to file 
        // N.B. - (this is the original BAM filename, method will modify it to use applicable extension)
        virtual bool Write(const std::string& bamFilename);
        
    // derived-classes MUST provide implementation
    protected:
        // clear all current index offset data in memory
        virtual void ClearAllData(void) =0;
        // return file position after header metadata
        virtual const off_t DataBeginOffset(void) const =0;
        // return true if all index data is cached
        virtual bool HasFullDataCache(void) const =0;
        // clears index data from all references except the first
        virtual void KeepOnlyFirstReferenceOffsets(void) =0;
        // load index data for all references, return true if loaded OK
        // @saveData - save data in memory if true, just read & discard if false
        virtual bool LoadAllReferences(bool saveData = true) =0;
        // load first reference from file, return true if loaded OK
        // @saveData - save data in memory if true, just read & discard if false
        virtual bool LoadFirstReference(bool saveData = true) =0;
        // load header data from index file, return true if loaded OK
        virtual bool LoadHeader(void) =0;
        // position file pointer to first reference begin, return true if skipped OK
        virtual bool SkipToFirstReference(void) =0;
        // write index reference data
        virtual bool WriteAllReferences(void) =0;
        // write index header data
        virtual bool WriteHeader(void) =0;

    // internal methods
    protected:
        // rewind index file to beginning of index data, return true if rewound OK
        bool Rewind(void);

    private:
        // return true if FILE* is open
        bool IsOpen(void) const;
        // opens index file according to requested mode, return true if opened OK
        bool OpenIndexFile(const std::string& filename, const std::string& mode);
        // updates in-memory cache of index data, depending on current cache mode
        void UpdateCache(void);

    // factory methods for returning proper BamIndex-derived type based on available index files
    public:
      
        // returns index based on BAM filename 'stub'
        // checks first for preferred type, returns that type if found
        // (if not found, attmempts to load other type(s), returns 0 if NONE found)
        //
        // ** default preferred type is BamToolsIndex ** use this anytime it exists
        enum PreferredIndexType { BAMTOOLS = 0, STANDARD };
        static BamIndex* FromBamFilename(const std::string&   bamFilename,
                                         BamTools::BgzfData*  bgzf,
                                         BamTools::BamReader* reader, 
                                         const BamIndex::PreferredIndexType& type = BamIndex::BAMTOOLS);
        
        // returns index based on explicitly named index file (or 0 if not found)
        static BamIndex* FromIndexFilename(const std::string&   indexFilename,
                                           BamTools::BgzfData*  bgzf,
                                           BamTools::BamReader* reader);

    // data members
    protected:
        BamTools::BgzfData*  m_BGZF;
        BamTools::BamReader* m_reader;
        BamTools::RefVector  m_references;
        BamIndex::BamIndexCacheMode m_cacheMode;
        FILE* m_indexStream;


    // friends
    friend class Internal::BamStandardIndex;
    friend class Internal::BamToolsIndex;
};

} // namespace BamTools

#endif // BAM_INDEX_H
