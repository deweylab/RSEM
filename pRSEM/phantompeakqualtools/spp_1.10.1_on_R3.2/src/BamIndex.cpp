// ***************************************************************************
// BamIndex.cpp (c) 2009 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 22 November 2010 (DB)
// ---------------------------------------------------------------------------
// Provides index functionality - both for the default (standardized) BAM 
// index format (.bai) as well as a BamTools-specific (nonstandard) index 
// format (.bti).
// ***************************************************************************

#include <BamIndex.h>
#include <BamReader.h>
#include <BGZF.h>
#include <BamStandardIndex_p.h>
#include <BamToolsIndex_p.h>
using namespace BamTools;
using namespace BamTools::Internal;

#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <map>
using namespace std;

// --------------------------------------------------
// BamIndex factory methods

// returns index based on BAM filename 'stub'
// checks first for preferred type, returns that type if found
// (if not found, attmempts to load other type(s), returns 0 if NONE found)
//
// ** default preferred type is BamToolsIndex ** use this anytime it exists
BamIndex* BamIndex::FromBamFilename(const std::string& bamFilename,
				    BamTools::BgzfData* bgzf,
				    BamTools::BamReader* reader,
				    const BamIndex::PreferredIndexType& type)
{
    // ---------------------------------------------------
    // attempt to load preferred type first

    const std::string bamtoolsIndexFilename = bamFilename + ".bti";
    const bool bamtoolsIndexExists = BamTools::FileExists(bamtoolsIndexFilename);
    if ( (type == BamIndex::BAMTOOLS) && bamtoolsIndexExists )
        return new BamToolsIndex(bgzf, reader);

    const std::string standardIndexFilename = bamFilename + ".bai";
    const bool standardIndexExists = BamTools::FileExists(standardIndexFilename);
    if ( (type == BamIndex::STANDARD) && standardIndexExists )
        return new BamStandardIndex(bgzf, reader);

    // ----------------------------------------------------
    // preferred type could not be found, try other (non-preferred) types
    // if none found, return 0

    if ( bamtoolsIndexExists ) return new BamToolsIndex(bgzf, reader);
    if ( standardIndexExists ) return new BamStandardIndex(bgzf, reader);
        return 0;
}

// returns index based on explicitly named index file (or 0 if not found)
BamIndex* BamIndex::FromIndexFilename(const std::string&   indexFilename,
				      BamTools::BgzfData*  bgzf,
				      BamTools::BamReader* reader)
{
    // see if specified file exists
    const bool indexExists = BamTools::FileExists(indexFilename);
    if ( !indexExists ) return 0;

    const std::string bamtoolsIndexExtension(".bti");
    const std::string standardIndexExtension(".bai");

    // if has bamtoolsIndexExtension
    if ( indexFilename.find(bamtoolsIndexExtension) == (indexFilename.length() - bamtoolsIndexExtension.length()) )
        return new BamToolsIndex(bgzf, reader);

     // if has standardIndexExtension
    if ( indexFilename.find(standardIndexExtension) == (indexFilename.length() - standardIndexExtension.length()) )
        return new BamStandardIndex(bgzf, reader);

    // otherwise, unsupported file type
    return 0;
}

// -------------------------------
// BamIndex implementation

// ctor
BamIndex::BamIndex(BamTools::BgzfData* bgzf, BamTools::BamReader* reader)
    : m_BGZF(bgzf)
    , m_reader(reader)
    , m_cacheMode(BamIndex::LimitedIndexCaching)
    , m_indexStream(0)
{ 
    if ( m_reader && m_reader->IsOpen() ) 
        m_references = m_reader->GetReferenceData();
}

// dtor
BamIndex::~BamIndex(void) {
    if ( IsOpen() )
        fclose(m_indexStream);
}

// return true if FILE* is open
bool BamIndex::IsOpen(void) const {
    return ( m_indexStream != 0 );
}

// loads existing data from file into memory
bool BamIndex::Load(const string& filename)  {

    // open index file, abort on error
    if ( !OpenIndexFile(filename, "rb") ) {
        fprintf(stderr, "ERROR: Unable to open the BAM index file %s for reading.\n", filename.c_str());
        return false;
    }

    // check magic number
    if ( !LoadHeader() ) {
        fclose(m_indexStream);
        return false;
    }

    // load reference data (but only keep in memory if full caching requested)
    bool saveInitialLoad = ( m_cacheMode == BamIndex::FullIndexCaching );
    if ( !LoadAllReferences(saveInitialLoad) ) {
        fclose(m_indexStream);
        return false;
    }

    // update index cache based on selected mode
    UpdateCache();

    // return success
    return true;
}

// opens index file for reading/writing, return true if opened OK
bool BamIndex::OpenIndexFile(const string& filename, const string& mode) {
    m_indexStream = fopen(filename.c_str(), mode.c_str());
    return ( m_indexStream != 0 );
}

// rewind index file to beginning of index data, return true if rewound OK
bool BamIndex::Rewind(void) {
    return ( fseek64(m_indexStream, DataBeginOffset(), SEEK_SET) == 0 );
}

// change the index caching behavior
void BamIndex::SetCacheMode(const BamIndexCacheMode mode) {
    if ( mode != m_cacheMode ) {
        m_cacheMode = mode;
        UpdateCache();
    }
}

// updates in-memory cache of index data, depending on current cache mode
void BamIndex::UpdateCache(void) {

    // skip if file not open
    if ( !IsOpen() ) return;

    // reflect requested cache mode behavior
    switch ( m_cacheMode ) {

        case (BamIndex::FullIndexCaching) :
            Rewind();
            LoadAllReferences(true);
            break;

        case (BamIndex::LimitedIndexCaching) :
            if ( HasFullDataCache() )
                KeepOnlyFirstReferenceOffsets();
            else {
                ClearAllData();
                SkipToFirstReference();
                LoadFirstReference(true);
            }
            break;
        case(BamIndex::NoIndexCaching) :
            ClearAllData();
            break;
        default :
            // unreachable
            ;
    }
}

// writes in-memory index data out to file
bool BamIndex::Write(const string& bamFilename) {

    // open index file for writing
    string indexFilename = bamFilename + Extension();
    if ( !OpenIndexFile(indexFilename, "wb") ) {
        fprintf(stderr, "ERROR: Could not open file to save index.\n");
        return false;
    }

    // write index header data
    if ( !WriteHeader() ) {
        fprintf(stderr, "ERROR: There was a problem writing index metadata to new index file.\n");
        fflush(m_indexStream);
        fclose(m_indexStream);
        exit(1);
    }

    // write main index data
    if ( !WriteAllReferences() ) {
        fprintf(stderr, "ERROR: There was a problem writing index data to new index file.\n");
        fflush(m_indexStream);
        fclose(m_indexStream);
        exit(1);
    }

    // flush any remaining output, rewind file, and return success
    fflush(m_indexStream);
    fclose(m_indexStream);

    // re-open index file for later reading
    if ( !OpenIndexFile(indexFilename, "rb") ) {
        fprintf(stderr, "ERROR: Could not open newly created index file for reading.\n");
        return false;
    }

    // return success/failure of write
    return true;
}
