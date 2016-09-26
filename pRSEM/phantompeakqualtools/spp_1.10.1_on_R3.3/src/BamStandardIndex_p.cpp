// ***************************************************************************
// BamStandardIndex.cpp (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 22 November 2010 (DB)
// ---------------------------------------------------------------------------
// Provides index operations for the standardized BAM index format (".bai")
// ***************************************************************************

#include <BamAlignment.h>
#include <BamReader.h>
#include <BGZF.h>
#include <BamStandardIndex_p.h>
using namespace BamTools;
using namespace BamTools::Internal;

#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <map>
using namespace std;

BamStandardIndex::BamStandardIndex(BgzfData* bgzf, BamReader* reader)
    : BamIndex(bgzf, reader)
    , m_dataBeginOffset(0)
    , m_hasFullDataCache(false)
{
    m_isBigEndian = BamTools::SystemIsBigEndian();
}

BamStandardIndex::~BamStandardIndex(void) {
    ClearAllData();
}

// calculate bins that overlap region
int BamStandardIndex::BinsFromRegion(const BamRegion& region,
				     const bool isRightBoundSpecified,
				     uint16_t bins[MAX_BIN])
{
    // get region boundaries
    uint32_t begin = (unsigned int)region.LeftPosition;
    uint32_t end;

    // if right bound specified AND left&right bounds are on same reference
    // OK to use right bound position
    if ( isRightBoundSpecified && ( region.LeftRefID == region.RightRefID ) )
	end = (unsigned int)region.RightPosition;

    // otherwise, use end of left bound reference as cutoff
    else
	end = (unsigned int)m_references.at(region.LeftRefID).RefLength - 1;

    // initialize list, bin '0' always a valid bin
    int i = 0;
    bins[i++] = 0;

    // get rest of bins that contain this region
    unsigned int k;
    for (k =    1 + (begin>>26); k <=    1 + (end>>26); ++k) { bins[i++] = k; }
    for (k =    9 + (begin>>23); k <=    9 + (end>>23); ++k) { bins[i++] = k; }
    for (k =   73 + (begin>>20); k <=   73 + (end>>20); ++k) { bins[i++] = k; }
    for (k =  585 + (begin>>17); k <=  585 + (end>>17); ++k) { bins[i++] = k; }
    for (k = 4681 + (begin>>14); k <= 4681 + (end>>14); ++k) { bins[i++] = k; }

    // return number of bins stored
    return i;
}

// creates index data (in-memory) from current reader data
bool BamStandardIndex::Build(void) {

    // be sure reader & BGZF file are valid & open for reading
    if ( m_reader == 0 || m_BGZF == 0 || !m_BGZF->IsOpen )
	return false;

    // move file pointer to beginning of alignments
    m_reader->Rewind();

    // get reference count, reserve index space
    const int numReferences = (int)m_references.size();
    m_indexData.clear();
    m_hasFullDataCache = false;
    SetReferenceCount(numReferences);

    // sets default constant for bin, ID, offset, coordinate variables
    const uint32_t defaultValue = 0xffffffffu;

    // bin data
    uint32_t saveBin(defaultValue);
    uint32_t lastBin(defaultValue);

    // reference ID data
    int32_t saveRefID(defaultValue);
    int32_t lastRefID(defaultValue);

    // offset data
    uint64_t saveOffset = m_BGZF->Tell();
    uint64_t lastOffset = saveOffset;

    // coordinate data
    int32_t lastCoordinate = defaultValue;

    BamAlignment bAlignment;
    while ( m_reader->GetNextAlignmentCore(bAlignment) ) {

	// change of chromosome, save ID, reset bin
	if ( lastRefID != bAlignment.RefID ) {
	    lastRefID = bAlignment.RefID;
	    lastBin   = defaultValue;
	}

	// if lastCoordinate greater than BAM position - file not sorted properly
	else if ( lastCoordinate > bAlignment.Position ) {
	    fprintf(stderr, "BAM file not properly sorted:\n");
	    fprintf(stderr, "Alignment %s : %d > %d on reference (id = %d)", bAlignment.Name.c_str(),
		    lastCoordinate, bAlignment.Position, bAlignment.RefID);
	    exit(1);
	}

	// if valid reference && BAM bin spans some minimum cutoff (smaller bin ids span larger regions)
	if ( (bAlignment.RefID >= 0) && (bAlignment.Bin < 4681) ) {

	    // save linear offset entry (matched to BAM entry refID)
	    BamStandardIndexData::iterator indexIter = m_indexData.find(bAlignment.RefID);
	    if ( indexIter == m_indexData.end() ) return false; // error
	    ReferenceIndex& refIndex = (*indexIter).second;
	    LinearOffsetVector& offsets = refIndex.Offsets;
	    SaveLinearOffset(offsets, bAlignment, lastOffset);
	}

	// if current BamAlignment bin != lastBin, "then possibly write the binning index"
	if ( bAlignment.Bin != lastBin ) {

	    // if not first time through
	    if ( saveBin != defaultValue ) {

		// save Bam bin entry
		BamStandardIndexData::iterator indexIter = m_indexData.find(saveRefID);
		if ( indexIter == m_indexData.end() ) return false; // error
		ReferenceIndex& refIndex = (*indexIter).second;
		BamBinMap& binMap = refIndex.Bins;
		SaveBinEntry(binMap, saveBin, saveOffset, lastOffset);
	    }

	    // update saveOffset
	    saveOffset = lastOffset;

	    // update bin values
	    saveBin = bAlignment.Bin;
	    lastBin = bAlignment.Bin;

	    // update saveRefID
	    saveRefID = bAlignment.RefID;

	    // if invalid RefID, break out
	    if ( saveRefID < 0 ) break;
	}

	// make sure that current file pointer is beyond lastOffset
	if ( m_BGZF->Tell() <= (int64_t)lastOffset ) {
	    fprintf(stderr, "Error in BGZF offsets.\n");
	    exit(1);
	}

	// update lastOffset
	lastOffset = m_BGZF->Tell();

	// update lastCoordinate
	lastCoordinate = bAlignment.Position;
    }

    // save any leftover BAM data (as long as refID is valid)
    if ( saveRefID >= 0 ) {
	// save Bam bin entry
	BamStandardIndexData::iterator indexIter = m_indexData.find(saveRefID);
	if ( indexIter == m_indexData.end() ) return false; // error
	ReferenceIndex& refIndex = (*indexIter).second;
	BamBinMap& binMap = refIndex.Bins;
	SaveBinEntry(binMap, saveBin, saveOffset, lastOffset);
    }

    // simplify index by merging chunks
    MergeChunks();

    // iterate through references in index
    // sort offsets in linear offset vector
    BamStandardIndexData::iterator indexIter = m_indexData.begin();
    BamStandardIndexData::iterator indexEnd  = m_indexData.end();
    for ( int i = 0; indexIter != indexEnd; ++indexIter, ++i ) {

	// get reference index data
	ReferenceIndex& refIndex = (*indexIter).second;
	LinearOffsetVector& offsets = refIndex.Offsets;

	// sort linear offsets
	sort(offsets.begin(), offsets.end());
    }

    // rewind file pointer to beginning of alignments, return success/fail
    return m_reader->Rewind();
}

// check index file magic number, return true if OK
bool BamStandardIndex::CheckMagicNumber(void) {

    // read in magic number
    char magic[4];
    size_t elementsRead = fread(magic, sizeof(char), 4, m_indexStream);

    // compare to expected value
    if ( strncmp(magic, "BAI\1", 4) != 0 ) {
	fprintf(stderr, "Problem with index file - invalid format.\n");
	fclose(m_indexStream);
	return false;
    }

    // return success/failure of load
    return (elementsRead == 4);
}

// clear all current index offset data in memory
void BamStandardIndex::ClearAllData(void) {
    BamStandardIndexData::const_iterator indexIter = m_indexData.begin();
    BamStandardIndexData::const_iterator indexEnd  = m_indexData.end();
    for ( ; indexIter != indexEnd; ++indexIter ) {
	const int& refId = (*indexIter).first;
	ClearReferenceOffsets(refId);
    }
}

// clear all index offset data for desired reference
void BamStandardIndex::ClearReferenceOffsets(const int& refId) {

    // look up refId, skip if not found
    BamStandardIndexData::iterator indexIter = m_indexData.find(refId);
    if ( indexIter == m_indexData.end() ) return ;

    // clear reference data
    ReferenceIndex& refEntry = (*indexIter).second;
    refEntry.Bins.clear();
    refEntry.Offsets.clear();

    // set flag
    m_hasFullDataCache = false;
}

// return file position after header metadata
const off_t BamStandardIndex::DataBeginOffset(void) const {
    return m_dataBeginOffset;
}

// calculates offset(s) for a given region
bool BamStandardIndex::GetOffsets(const BamRegion& region,
				  const bool isRightBoundSpecified,
				  vector<int64_t>& offsets,
				  bool* hasAlignmentsInRegion)
{
    // return false if leftBound refID is not found in index data
    if ( m_indexData.find(region.LeftRefID) == m_indexData.end() )
	return false;

    // load index data for region if not already cached
    if ( !IsDataLoaded(region.LeftRefID) ) {
	bool loadedOk = true;
	loadedOk &= SkipToReference(region.LeftRefID);
	loadedOk &= LoadReference(region.LeftRefID);
	if ( !loadedOk ) return false;
    }

    // calculate which bins overlap this region
    uint16_t* bins = (uint16_t*)calloc(MAX_BIN, 2);
    int numBins = BinsFromRegion(region, isRightBoundSpecified, bins);

    // get bins for this reference
    BamStandardIndexData::const_iterator indexIter = m_indexData.find(region.LeftRefID);
    if ( indexIter == m_indexData.end() ) return false; // error
    const ReferenceIndex& refIndex = (*indexIter).second;
    const BamBinMap& binMap        = refIndex.Bins;

    // get minimum offset to consider
    const LinearOffsetVector& linearOffsets = refIndex.Offsets;
    const uint64_t minOffset = ( (unsigned int)(region.LeftPosition>>BAM_LIDX_SHIFT) >= linearOffsets.size() )
			       ? 0 : linearOffsets.at(region.LeftPosition>>BAM_LIDX_SHIFT);

    // store all alignment 'chunk' starts (file offsets) for bins in this region
    for ( int i = 0; i < numBins; ++i ) {

	const uint16_t binKey = bins[i];
	map<uint32_t, ChunkVector>::const_iterator binIter = binMap.find(binKey);
	if ( (binIter != binMap.end()) && ((*binIter).first == binKey) ) {

	    // iterate over chunks
	    const ChunkVector& chunks = (*binIter).second;
	    std::vector<Chunk>::const_iterator chunksIter = chunks.begin();
	    std::vector<Chunk>::const_iterator chunksEnd  = chunks.end();
	    for ( ; chunksIter != chunksEnd; ++chunksIter) {

		// if valid chunk found, store its file offset
		const Chunk& chunk = (*chunksIter);
		if ( chunk.Stop > minOffset )
		    offsets.push_back( chunk.Start );
	    }
	}
    }

    // clean up memory
    free(bins);

    // sort the offsets before returning
    sort(offsets.begin(), offsets.end());

    // set flag & return success
    *hasAlignmentsInRegion = (offsets.size() != 0 );

    // if cache mode set to none, dump the data we just loaded
    if (m_cacheMode == BamIndex::NoIndexCaching )
	ClearReferenceOffsets(region.LeftRefID);

    // return succes
    return true;
}

// returns whether reference has alignments or no
bool BamStandardIndex::HasAlignments(const int& refId) const {
    BamStandardIndexData::const_iterator indexIter = m_indexData.find(refId);
    if ( indexIter == m_indexData.end() ) return false; // error
    const ReferenceIndex& refEntry = (*indexIter).second;
    return refEntry.HasAlignments;
}

// return true if all index data is cached
bool BamStandardIndex::HasFullDataCache(void) const {
    return m_hasFullDataCache;
}

// returns true if index cache has data for desired reference
bool BamStandardIndex::IsDataLoaded(const int& refId) const {

    // look up refId, return false if not found
    BamStandardIndexData::const_iterator indexIter = m_indexData.find(refId);
    if ( indexIter == m_indexData.end() ) return false;

    // see if reference has alignments
    // if not, it's not a problem to have no offset data
    const ReferenceIndex& refEntry = (*indexIter).second;
    if ( !refEntry.HasAlignments ) return true;

    // return whether bin map contains data
    return ( !refEntry.Bins.empty() );
}

// attempts to use index to jump to region; returns success/fail
bool BamStandardIndex::Jump(const BamRegion& region, bool* hasAlignmentsInRegion) {

    // be sure reader & BGZF file are valid & open for reading
    if ( m_reader == 0 || m_BGZF == 0 || !m_BGZF->IsOpen )
	return false;

    // make sure left-bound position is valid
    if ( region.LeftPosition > m_references.at(region.LeftRefID).RefLength )
	return false;

    // calculate offsets for this region
    // if failed, print message, set flag, and return failure
    vector<int64_t> offsets;
    if ( !GetOffsets(region, region.isRightBoundSpecified(), offsets, hasAlignmentsInRegion) ) {
	fprintf(stderr, "ERROR: Could not jump: unable to calculate offset(s) for specified region.\n");
	*hasAlignmentsInRegion = false;
	return false;
    }

    // iterate through offsets
    BamAlignment bAlignment;
    bool result = true;
    for ( vector<int64_t>::const_iterator o = offsets.begin(); o != offsets.end(); ++o) {

	// attempt seek & load first available alignment
	// set flag to true if data exists
	result &= m_BGZF->Seek(*o);
	*hasAlignmentsInRegion = m_reader->GetNextAlignmentCore(bAlignment);

	// if this alignment corresponds to desired position
	// return success of seeking back to the offset before the 'current offset' (to cover overlaps)
	if ( ((bAlignment.RefID == region.LeftRefID) &&
	       ((bAlignment.Position + bAlignment.Length) > region.LeftPosition)) ||
	     (bAlignment.RefID > region.LeftRefID) )
	{
	    if ( o != offsets.begin() ) --o;
	    return m_BGZF->Seek(*o);
	}
    }

    // if error in jumping, print message & set flag
    if ( !result ) {
	fprintf(stderr, "ERROR: Could not jump: unable to determine correct offset for specified region.\n");
	*hasAlignmentsInRegion = false;
    }

    // return success/failure
    return result;
}

// clears index data from all references except the first
void BamStandardIndex::KeepOnlyFirstReferenceOffsets(void) {
    BamStandardIndexData::const_iterator indexBegin = m_indexData.begin();
    KeepOnlyReferenceOffsets((*indexBegin).first);
}

// clears index data from all references except the one specified
void BamStandardIndex::KeepOnlyReferenceOffsets(const int& refId) {
    BamStandardIndexData::iterator mapIter = m_indexData.begin();
    BamStandardIndexData::iterator mapEnd  = m_indexData.end();
    for ( ; mapIter != mapEnd; ++mapIter ) {
	const int entryRefId = (*mapIter).first;
	if ( entryRefId != refId )
	    ClearReferenceOffsets(entryRefId);
    }
}

bool BamStandardIndex::LoadAllReferences(bool saveData) {

    // skip if data already loaded
    if ( m_hasFullDataCache ) return true;

    // get number of reference sequences
    uint32_t numReferences;
    if ( !LoadReferenceCount((int&)numReferences) )
	return false;

    // iterate over reference entries
    bool loadedOk = true;
    for ( int i = 0; i < (int)numReferences; ++i )
	loadedOk &= LoadReference(i, saveData);

    // set flag
    if ( loadedOk && saveData )
	m_hasFullDataCache = true;

    // return success/failure of loading references
    return loadedOk;
}

// load header data from index file, return true if loaded OK
bool BamStandardIndex::LoadHeader(void) {

    bool loadedOk = CheckMagicNumber();

    // store offset of beginning of data
    m_dataBeginOffset = ftell64(m_indexStream);

    // return success/failure of load
    return loadedOk;
}

// load a single index bin entry from file, return true if loaded OK
// @saveData - save data in memory if true, just read & discard if false
bool BamStandardIndex::LoadBin(ReferenceIndex& refEntry, bool saveData) {

    size_t elementsRead = 0;

    // get bin ID
    uint32_t binId;
    elementsRead += fread(&binId, sizeof(binId), 1, m_indexStream);
    if ( m_isBigEndian ) SwapEndian_32(binId);

    // load alignment chunks for this bin
    ChunkVector chunks;
    bool chunksOk = LoadChunks(chunks, saveData);

    // store bin entry
    if ( chunksOk && saveData )
	refEntry.Bins.insert(pair<uint32_t, ChunkVector>(binId, chunks));

    // return success/failure of load
    return ( (elementsRead == 1) && chunksOk );
}

bool BamStandardIndex::LoadBins(ReferenceIndex& refEntry, bool saveData) {

    size_t elementsRead = 0;

    // get number of bins
    int32_t numBins;
    elementsRead += fread(&numBins, sizeof(numBins), 1, m_indexStream);
    if ( m_isBigEndian ) SwapEndian_32(numBins);

    // set flag
    refEntry.HasAlignments = ( numBins != 0 );

    // iterate over bins
    bool binsOk = true;
    for ( int i = 0; i < numBins; ++i )
	binsOk &= LoadBin(refEntry, saveData);

    // return success/failure of load
    return ( (elementsRead == 1) && binsOk );
}

// load a single index bin entry from file, return true if loaded OK
// @saveData - save data in memory if true, just read & discard if false
bool BamStandardIndex::LoadChunk(ChunkVector& chunks, bool saveData) {

    size_t elementsRead = 0;

    // read in chunk data
    uint64_t start;
    uint64_t stop;
    elementsRead += fread(&start, sizeof(start), 1, m_indexStream);
    elementsRead += fread(&stop,  sizeof(stop),  1, m_indexStream);

    // swap endian-ness if necessary
    if ( m_isBigEndian ) {
	SwapEndian_64(start);
	SwapEndian_64(stop);
    }

    // save data if requested
    if ( saveData ) chunks.push_back( Chunk(start, stop) );

    // return success/failure of load
    return ( elementsRead == 2 );
}

bool BamStandardIndex::LoadChunks(ChunkVector& chunks, bool saveData) {

    size_t elementsRead = 0;

    // read in number of chunks
    uint32_t numChunks;
    elementsRead += fread(&numChunks, sizeof(numChunks), 1, m_indexStream);
    if ( m_isBigEndian ) SwapEndian_32(numChunks);

    // initialize space for chunks if we're storing this data
    if ( saveData ) chunks.reserve(numChunks);

    // iterate over chunks
    bool chunksOk = true;
    for ( int i = 0; i < (int)numChunks; ++i )
	chunksOk &= LoadChunk(chunks, saveData);

    // sort chunk vector
    sort( chunks.begin(), chunks.end(), ChunkLessThan );

    // return success/failure of load
    return ( (elementsRead == 1) && chunksOk );
}

// load a single index linear offset entry from file, return true if loaded OK
// @saveData - save data in memory if true, just read & discard if false
bool BamStandardIndex::LoadLinearOffsets(ReferenceIndex& refEntry, bool saveData) {

    size_t elementsRead = 0;

    // read in number of linear offsets
    int32_t numLinearOffsets;
    elementsRead += fread(&numLinearOffsets, sizeof(numLinearOffsets), 1, m_indexStream);
    if ( m_isBigEndian ) SwapEndian_32(numLinearOffsets);

    // set up destination vector (if we're saving the data)
    LinearOffsetVector linearOffsets;
    if ( saveData ) linearOffsets.reserve(numLinearOffsets);

    // iterate over linear offsets
    uint64_t linearOffset;
    for ( int i = 0; i < numLinearOffsets; ++i ) {
	elementsRead += fread(&linearOffset, sizeof(linearOffset), 1, m_indexStream);
	if ( m_isBigEndian ) SwapEndian_64(linearOffset);
	if ( saveData ) linearOffsets.push_back(linearOffset);
    }

    // sort linear offsets
    sort ( linearOffsets.begin(), linearOffsets.end() );

    // save in reference index entry if desired
    if ( saveData ) refEntry.Offsets = linearOffsets;

    // return success/failure of load
    return ( elementsRead == (size_t)(numLinearOffsets + 1) );
}

bool BamStandardIndex::LoadFirstReference(bool saveData) {
    BamStandardIndexData::const_iterator indexBegin = m_indexData.begin();
    return LoadReference((*indexBegin).first, saveData);
}

// load a single reference from file, return true if loaded OK
// @saveData - save data in memory if true, just read & discard if false
bool BamStandardIndex::LoadReference(const int& refId, bool saveData) {

    // look up refId
    BamStandardIndexData::iterator indexIter = m_indexData.find(refId);

    // if reference not previously loaded, create new entry
    if ( indexIter == m_indexData.end() ) {
	ReferenceIndex newEntry;
	newEntry.HasAlignments = false;
	m_indexData.insert( pair<int32_t, ReferenceIndex>(refId, newEntry) );
    }

    // load reference data
    indexIter = m_indexData.find(refId);
    ReferenceIndex& entry = (*indexIter).second;
    bool loadedOk = true;
    loadedOk &= LoadBins(entry, saveData);
    loadedOk &= LoadLinearOffsets(entry, saveData);
    return loadedOk;
}

// loads number of references, return true if loaded OK
bool BamStandardIndex::LoadReferenceCount(int& numReferences) {

    size_t elementsRead = 0;

    // read reference count
    elementsRead += fread(&numReferences, sizeof(numReferences), 1, m_indexStream);
    if ( m_isBigEndian ) SwapEndian_32(numReferences);

    // return success/failure of load
    return ( elementsRead == 1 );
}

// merges 'alignment chunks' in BAM bin (used for index building)
void BamStandardIndex::MergeChunks(void) {

    // iterate over reference enties
    BamStandardIndexData::iterator indexIter = m_indexData.begin();
    BamStandardIndexData::iterator indexEnd  = m_indexData.end();
    for ( ; indexIter != indexEnd; ++indexIter ) {

	// get BAM bin map for this reference
	ReferenceIndex& refIndex = (*indexIter).second;
	BamBinMap& bamBinMap = refIndex.Bins;

	// iterate over BAM bins
	BamBinMap::iterator binIter = bamBinMap.begin();
	BamBinMap::iterator binEnd  = bamBinMap.end();
	for ( ; binIter != binEnd; ++binIter ) {

	    // get chunk vector for this bin
	    ChunkVector& binChunks = (*binIter).second;
	    if ( binChunks.size() == 0 ) continue;

	    ChunkVector mergedChunks;
	    mergedChunks.push_back( binChunks[0] );

	    // iterate over chunks
	    int i = 0;
	    ChunkVector::iterator chunkIter = binChunks.begin();
	    ChunkVector::iterator chunkEnd  = binChunks.end();
	    for ( ++chunkIter; chunkIter != chunkEnd; ++chunkIter) {

		// get 'currentChunk' based on numeric index
		Chunk& currentChunk = mergedChunks[i];

		// get iteratorChunk based on vector iterator
		Chunk& iteratorChunk = (*chunkIter);

		// if chunk ends where (iterator) chunk starts, then merge
		if ( currentChunk.Stop>>16 == iteratorChunk.Start>>16 )
		    currentChunk.Stop = iteratorChunk.Stop;

		// otherwise
		else {
		    // set currentChunk + 1 to iteratorChunk
		    mergedChunks.push_back(iteratorChunk);
		    ++i;
		}
	    }

	    // saved merged chunk vector
	    (*binIter).second = mergedChunks;
	}
    }
}

// saves BAM bin entry for index
void BamStandardIndex::SaveBinEntry(BamBinMap& binMap,
				    const uint32_t& saveBin,
				    const uint64_t& saveOffset,
				    const uint64_t& lastOffset)
{
    // look up saveBin
    BamBinMap::iterator binIter = binMap.find(saveBin);

    // create new chunk
    Chunk newChunk(saveOffset, lastOffset);

    // if entry doesn't exist
    if ( binIter == binMap.end() ) {
	ChunkVector newChunks;
	newChunks.push_back(newChunk);
	binMap.insert( pair<uint32_t, ChunkVector>(saveBin, newChunks));
    }

    // otherwise
    else {
	ChunkVector& binChunks = (*binIter).second;
	binChunks.push_back( newChunk );
    }
}

// saves linear offset entry for index
void BamStandardIndex::SaveLinearOffset(LinearOffsetVector& offsets,
					const BamAlignment& bAlignment,
					const uint64_t&     lastOffset)
{
    // get converted offsets
    int beginOffset = bAlignment.Position >> BAM_LIDX_SHIFT;
    int endOffset   = (bAlignment.GetEndPosition() - 1) >> BAM_LIDX_SHIFT;

    // resize vector if necessary
    int oldSize = offsets.size();
    int newSize = endOffset + 1;
    if ( oldSize < newSize )
	offsets.resize(newSize, 0);

    // store offset
    for( int i = beginOffset + 1; i <= endOffset; ++i ) {
	if ( offsets[i] == 0 )
	    offsets[i] = lastOffset;
    }
}

// initializes index data structure to hold @count references
void BamStandardIndex::SetReferenceCount(const int& count) {
    for ( int i = 0; i < count; ++i )
	m_indexData[i].HasAlignments = false;
}

bool BamStandardIndex::SkipToFirstReference(void) {
    BamStandardIndexData::const_iterator indexBegin = m_indexData.begin();
    return SkipToReference( (*indexBegin).first );
}

// position file pointer to desired reference begin, return true if skipped OK
bool BamStandardIndex::SkipToReference(const int& refId) {

    // attempt rewind
    if ( !Rewind() ) return false;

    // read in number of references
    uint32_t numReferences;
    size_t elementsRead = fread(&numReferences, sizeof(numReferences), 1, m_indexStream);
    if ( elementsRead != 1 ) return false;
    if ( m_isBigEndian ) SwapEndian_32(numReferences);

    // iterate over reference entries
    bool skippedOk = true;
    int currentRefId = 0;
    while (currentRefId != refId) {
	skippedOk &= LoadReference(currentRefId, false);
	++currentRefId;
    }

    // return success
    return skippedOk;
}

// write header to new index file
bool BamStandardIndex::WriteHeader(void) {

    size_t elementsWritten = 0;

    // write magic number
    elementsWritten += fwrite("BAI\1", sizeof(char), 4, m_indexStream);

    // store offset of beginning of data
    m_dataBeginOffset = ftell64(m_indexStream);

    // return success/failure of write
    return (elementsWritten == 4);
}

// write index data for all references to new index file
bool BamStandardIndex::WriteAllReferences(void) {

    size_t elementsWritten = 0;

    // write number of reference sequences
    int32_t numReferenceSeqs = m_indexData.size();
    if ( m_isBigEndian ) SwapEndian_32(numReferenceSeqs);
    elementsWritten += fwrite(&numReferenceSeqs, sizeof(numReferenceSeqs), 1, m_indexStream);

    // iterate over reference sequences
    bool refsOk = true;
    BamStandardIndexData::const_iterator indexIter = m_indexData.begin();
    BamStandardIndexData::const_iterator indexEnd  = m_indexData.end();
    for ( ; indexIter != indexEnd; ++ indexIter )
	refsOk &= WriteReference( (*indexIter).second );

    // return success/failure of write
    return ( (elementsWritten == 1) && refsOk );
}

// write index data for bin to new index file
bool BamStandardIndex::WriteBin(const uint32_t& binId, const ChunkVector& chunks) {

    size_t elementsWritten = 0;

    // write BAM bin ID
    uint32_t binKey = binId;
    if ( m_isBigEndian ) SwapEndian_32(binKey);
    elementsWritten += fwrite(&binKey, sizeof(binKey), 1, m_indexStream);

    // write chunks
    bool chunksOk = WriteChunks(chunks);

    // return success/failure of write
    return ( (elementsWritten == 1) && chunksOk );
}

// write index data for bins to new index file
bool BamStandardIndex::WriteBins(const BamBinMap& bins) {

    size_t elementsWritten = 0;

    // write number of bins
    int32_t binCount = bins.size();
    if ( m_isBigEndian ) SwapEndian_32(binCount);
    elementsWritten += fwrite(&binCount, sizeof(binCount), 1, m_indexStream);

    // iterate over bins
    bool binsOk = true;
    BamBinMap::const_iterator binIter = bins.begin();
    BamBinMap::const_iterator binEnd  = bins.end();
    for ( ; binIter != binEnd; ++binIter )
	binsOk &= WriteBin( (*binIter).first, (*binIter).second );

    // return success/failure of write
    return ( (elementsWritten == 1) && binsOk );
}

// write index data for chunk entry to new index file
bool BamStandardIndex::WriteChunk(const Chunk& chunk) {

    size_t elementsWritten = 0;

    // localize alignment chunk offsets
    uint64_t start = chunk.Start;
    uint64_t stop  = chunk.Stop;

    // swap endian-ness if necessary
    if ( m_isBigEndian ) {
	SwapEndian_64(start);
	SwapEndian_64(stop);
    }

    // write to index file
    elementsWritten += fwrite(&start, sizeof(start), 1, m_indexStream);
    elementsWritten += fwrite(&stop,  sizeof(stop),  1, m_indexStream);

    // return success/failure of write
    return ( elementsWritten == 2 );
}

// write index data for chunk entry to new index file
bool BamStandardIndex::WriteChunks(const ChunkVector& chunks) {

    size_t elementsWritten = 0;

    // write chunks
    int32_t chunkCount = chunks.size();
    if ( m_isBigEndian ) SwapEndian_32(chunkCount);
    elementsWritten += fwrite(&chunkCount, sizeof(chunkCount), 1, m_indexStream);

    // iterate over chunks
    bool chunksOk = true;
    ChunkVector::const_iterator chunkIter = chunks.begin();
    ChunkVector::const_iterator chunkEnd  = chunks.end();
    for ( ; chunkIter != chunkEnd; ++chunkIter )
	chunksOk &= WriteChunk( (*chunkIter) );

    // return success/failure of write
    return ( (elementsWritten == 1) && chunksOk );
}

// write index data for linear offsets entry to new index file
bool BamStandardIndex::WriteLinearOffsets(const LinearOffsetVector& offsets) {

    size_t elementsWritten = 0;

    // write number of linear offsets
    int32_t offsetCount = offsets.size();
    if ( m_isBigEndian ) SwapEndian_32(offsetCount);
    elementsWritten += fwrite(&offsetCount, sizeof(offsetCount), 1, m_indexStream);

    // iterate over linear offsets
    LinearOffsetVector::const_iterator offsetIter = offsets.begin();
    LinearOffsetVector::const_iterator offsetEnd  = offsets.end();
    for ( ; offsetIter != offsetEnd; ++offsetIter ) {

	// write linear offset
	uint64_t linearOffset = (*offsetIter);
	if ( m_isBigEndian ) SwapEndian_64(linearOffset);
	elementsWritten += fwrite(&linearOffset, sizeof(linearOffset), 1, m_indexStream);
    }

    // return success/failure of write
    return ( elementsWritten == (size_t)(offsetCount + 1) );
}

// write index data for a single reference to new index file
bool BamStandardIndex::WriteReference(const ReferenceIndex& refEntry) {
    bool refOk = true;
    refOk &= WriteBins(refEntry.Bins);
    refOk &= WriteLinearOffsets(refEntry.Offsets);
    return refOk;
}
