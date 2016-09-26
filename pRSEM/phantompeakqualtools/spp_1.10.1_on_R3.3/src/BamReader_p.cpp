// ***************************************************************************
// BamReader_p.cpp (c) 2009 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 22 November 2010 (DB)
// ---------------------------------------------------------------------------
// Provides the basic functionality for reading BAM files
// ***************************************************************************

#include <BamReader.h>
#include <BGZF.h>
#include <BamReader_p.h>
#include <BamStandardIndex_p.h>
#include <BamToolsIndex_p.h>
using namespace BamTools;
using namespace BamTools::Internal;

#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>
using namespace std;

// constructor
BamReaderPrivate::BamReaderPrivate(BamReader* parent)
    : HeaderText("")
    , Index(0)
    , HasIndex(false)
    , AlignmentsBeginOffset(0)
//    , m_header(0)
    , IndexCacheMode(BamIndex::LimitedIndexCaching)
    , HasAlignmentsInRegion(true)
    , Parent(parent)
    , DNA_LOOKUP("=ACMGRSVTWYHKDBN")
    , CIGAR_LOOKUP("MIDNSHP")
{
    IsBigEndian = SystemIsBigEndian();
}

// destructor
BamReaderPrivate::~BamReaderPrivate(void) {
    Close();
}

// adjusts requested region if necessary (depending on where data actually begins)
void BamReaderPrivate::AdjustRegion(BamRegion& region) {

    // check for valid index first
    if ( Index == 0 ) return;

    // see if any references in region have alignments
    HasAlignmentsInRegion = false;
    int currentId = region.LeftRefID;

    const int rightBoundRefId = ( region.isRightBoundSpecified() ? region.RightRefID : References.size() - 1 );
    while ( currentId <= rightBoundRefId ) {
	HasAlignmentsInRegion = Index->HasAlignments(currentId);
	if ( HasAlignmentsInRegion ) break;
	++currentId;
    }

    // if no data found on any reference in region
    if ( !HasAlignmentsInRegion ) return;

    // if left bound of desired region had no data, use first reference that had data
    // otherwise, leave requested region as-is
    if ( currentId != region.LeftRefID ) {
	region.LeftRefID = currentId;
	region.LeftPosition = 0;
    }
}

// fills out character data for BamAlignment data
bool BamReaderPrivate::BuildCharData(BamAlignment& bAlignment) {

    // calculate character lengths/offsets
    const unsigned int dataLength     = bAlignment.SupportData.BlockLength - BAM_CORE_SIZE;
    const unsigned int seqDataOffset  = bAlignment.SupportData.QueryNameLength + (bAlignment.SupportData.NumCigarOperations * 4);
    const unsigned int qualDataOffset = seqDataOffset + (bAlignment.SupportData.QuerySequenceLength+1)/2;
    const unsigned int tagDataOffset  = qualDataOffset + bAlignment.SupportData.QuerySequenceLength;
    const unsigned int tagDataLength  = dataLength - tagDataOffset;

    // check offsets to see what char data exists
    const bool hasSeqData  = ( seqDataOffset  < dataLength );
    const bool hasQualData = ( qualDataOffset < dataLength );
    const bool hasTagData  = ( tagDataOffset  < dataLength );

    // set up char buffers
    const char* allCharData = bAlignment.SupportData.AllCharData.data();
    const char* seqData     = ( hasSeqData  ? (((const char*)allCharData) + seqDataOffset)  : (const char*)0 );
    const char* qualData    = ( hasQualData ? (((const char*)allCharData) + qualDataOffset) : (const char*)0 );
          char* tagData     = ( hasTagData  ? (((char*)allCharData) + tagDataOffset)        : (char*)0 );

    // store alignment name (relies on null char in name as terminator)
    bAlignment.Name.assign((const char*)(allCharData));

    // save query sequence
    bAlignment.QueryBases.clear();
    if ( hasSeqData ) {
        bAlignment.QueryBases.reserve(bAlignment.SupportData.QuerySequenceLength);
        for (unsigned int i = 0; i < bAlignment.SupportData.QuerySequenceLength; ++i) {
            char singleBase = DNA_LOOKUP[ ( (seqData[(i/2)] >> (4*(1-(i%2)))) & 0xf ) ];
            bAlignment.QueryBases.append(1, singleBase);
        }
    }

    // save qualities, converting from numeric QV to 'FASTQ-style' ASCII character
    bAlignment.Qualities.clear();
    if ( hasQualData ) {
        bAlignment.Qualities.reserve(bAlignment.SupportData.QuerySequenceLength);
        for (unsigned int i = 0; i < bAlignment.SupportData.QuerySequenceLength; ++i) {
            char singleQuality = (char)(qualData[i]+33);
            bAlignment.Qualities.append(1, singleQuality);
        }
    }

    // if QueryBases is empty (and this is a allowed case)
    if ( bAlignment.QueryBases.empty() )
        bAlignment.AlignedBases = bAlignment.QueryBases;

    // if QueryBases contains data, then build AlignedBases using CIGAR data
    else {

        // resize AlignedBases
        bAlignment.AlignedBases.clear();
        bAlignment.AlignedBases.reserve(bAlignment.SupportData.QuerySequenceLength);

        // iterate over CigarOps
        int k = 0;
        vector<CigarOp>::const_iterator cigarIter = bAlignment.CigarData.begin();
        vector<CigarOp>::const_iterator cigarEnd  = bAlignment.CigarData.end();
        for ( ; cigarIter != cigarEnd; ++cigarIter ) {

            const CigarOp& op = (*cigarIter);
            switch(op.Type) {

            case ('M') :
            case ('I') :
                bAlignment.AlignedBases.append(bAlignment.QueryBases.substr(k, op.Length)); // for 'M', 'I' - write bases
                // fall through

            case ('S') :
                k += op.Length;                                     // for 'S' - soft clip, skip over query bases
                break;

            case ('D') :
                bAlignment.AlignedBases.append(op.Length, '-');     // for 'D' - write gap character
                break;

            case ('P') :
                bAlignment.AlignedBases.append( op.Length, '*' );   // for 'P' - write padding character
                break;

            case ('N') :
                bAlignment.AlignedBases.append( op.Length, 'N' );  // for 'N' - write N's, skip bases in original query sequence
                break;

            case ('H') :
                break;  // for 'H' - hard clip, do nothing to AlignedBases, move to next op

            default:
                fprintf(stderr, "ERROR: Invalid Cigar op type\n"); // shouldn't get here
                exit(1);
            }
        }
    }

    // save tag data
    bAlignment.TagData.clear();
    if ( hasTagData ) {
        if ( IsBigEndian ) {
            int i = 0;
            while ( (unsigned int)i < tagDataLength ) {

                i += 2; // skip tag type (e.g. "RG", "NM", etc)
                uint8_t type = toupper(tagData[i]);     // lower & upper case letters have same meaning
                ++i;                                    // skip value type

                switch (type) {

                    case('A') :
                    case('C') :
                        ++i;
                        break;

                    case('S') :
                        SwapEndian_16p(&tagData[i]);
                        i += sizeof(uint16_t);
                        break;

                    case('F') :
                    case('I') :
                        SwapEndian_32p(&tagData[i]);
                        i += sizeof(uint32_t);
                        break;

                    case('D') :
                        SwapEndian_64p(&tagData[i]);
                        i += sizeof(uint64_t);
                        break;

                    case('H') :
                    case('Z') :
                        while (tagData[i]) { ++i; }
                        ++i; // increment one more for null terminator
                        break;

                    default :
                        fprintf(stderr, "ERROR: Invalid tag value type\n"); // shouldn't get here
                        exit(1);
                }
            }
        }

        // store tagData in alignment
        bAlignment.TagData.resize(tagDataLength);
        memcpy((char*)bAlignment.TagData.data(), tagData, tagDataLength);
    }

    // clear the core-only flag
    bAlignment.SupportData.HasCoreOnly = false;

    // return success
    return true;
}

// clear index data structure
void BamReaderPrivate::ClearIndex(void) {
    delete Index;
    Index = 0;
    HasIndex = false;
}

// closes the BAM file
void BamReaderPrivate::Close(void) {

    // close BGZF file stream
    mBGZF.Close();

    // clear out index data
    ClearIndex();

    // clear out header data
    HeaderText.clear();
//    if ( m_header ) {
//	delete m_header;
//	m_header = 0;
//    }

    // clear out region flags
    Region.clear();
}

// creates index for BAM file, saves to file
// default behavior is to create the BAM standard index (".bai")
// set flag to false to create the BamTools-specific index (".bti")
bool BamReaderPrivate::CreateIndex(bool useStandardIndex) {

    // clear out prior index data
    ClearIndex();

    // create index based on type requested
    if ( useStandardIndex )
	Index = new BamStandardIndex(&mBGZF, Parent);
    else
	Index = new BamToolsIndex(&mBGZF, Parent);

    // set index cache mode to full for writing
    Index->SetCacheMode(BamIndex::FullIndexCaching);

    // build new index
    bool ok = true;
    ok &= Index->Build();
    HasIndex = ok;

    // mark empty references
    MarkReferences();

    // attempt to save index data to file
    ok &= Index->Write(Filename);

    // set client's desired index cache mode
    Index->SetCacheMode(IndexCacheMode);

    // return success/fail of both building & writing index
    return ok;
}

const string BamReaderPrivate::GetHeaderText(void) const {

    return HeaderText;

//    if ( m_header )
//	return m_header->Text();
//    else
//	return string("");
}

// get next alignment (from specified region, if given)
bool BamReaderPrivate::GetNextAlignment(BamAlignment& bAlignment) {

    // if valid alignment found, attempt to parse char data, and return success/failure
    if ( GetNextAlignmentCore(bAlignment) )
	return BuildCharData(bAlignment);

    // no valid alignment found
    else return false;
}

// retrieves next available alignment core data (returns success/fail)
// ** DOES NOT parse any character data (read name, bases, qualities, tag data)
//    these can be accessed, if necessary, from the supportData
// useful for operations requiring ONLY positional or other alignment-related information
bool BamReaderPrivate::GetNextAlignmentCore(BamAlignment& bAlignment) {

    // if region is set but has no alignments
    if ( !Region.isNull() && !HasAlignmentsInRegion )
	return false;

    // if valid alignment available
    if ( LoadNextAlignment(bAlignment) ) {

	// set core-only flag
	bAlignment.SupportData.HasCoreOnly = true;

	// if region not specified with at least a left boundary, return success
	if ( !Region.isLeftBoundSpecified() ) return true;

	// determine region state (before, within, after)
	BamReaderPrivate::RegionState state = IsOverlap(bAlignment);

	// if alignment lies after region, return false
	if ( state == AFTER_REGION ) return false;

	while ( state != WITHIN_REGION ) {
	    // if no valid alignment available (likely EOF) return failure
	    if ( !LoadNextAlignment(bAlignment) ) return false;
	    // if alignment lies after region, return false (no available read within region)
	    state = IsOverlap(bAlignment);
	    if ( state == AFTER_REGION ) return false;
	}

	// return success (alignment found that overlaps region)
	return true;
    }

    // no valid alignment
    else return false;
}

// returns RefID for given RefName (returns References.size() if not found)
int BamReaderPrivate::GetReferenceID(const string& refName) const {

    // retrieve names from reference data
    vector<string> refNames;
    RefVector::const_iterator refIter = References.begin();
    RefVector::const_iterator refEnd  = References.end();
    for ( ; refIter != refEnd; ++refIter)
	refNames.push_back( (*refIter).RefName );

    // return 'index-of' refName ( if not found, returns refNames.size() )
    return distance(refNames.begin(), find(refNames.begin(), refNames.end(), refName));
}

// returns region state - whether alignment ends before, overlaps, or starts after currently specified region
// this *internal* method should ONLY called when (at least) IsLeftBoundSpecified == true
BamReaderPrivate::RegionState BamReaderPrivate::IsOverlap(BamAlignment& bAlignment) {

    // if alignment is on any reference sequence before left bound
    if ( bAlignment.RefID < Region.LeftRefID ) return BEFORE_REGION;

    // if alignment starts on left bound reference
    else if ( bAlignment.RefID == Region.LeftRefID ) {

	// if alignment starts at or after left boundary
	if ( bAlignment.Position >= Region.LeftPosition) {

	    // if right boundary is specified AND
	    // left/right boundaries are on same reference AND
	    // alignment starts past right boundary
	    if ( Region.isRightBoundSpecified() &&
		 Region.LeftRefID == Region.RightRefID &&
		 bAlignment.Position > Region.RightPosition )
		return AFTER_REGION;

	    // otherwise, alignment is within region
	    return WITHIN_REGION;
	}

	// alignment starts before left boundary
	else {
	    // check if alignment overlaps left boundary
	    if ( bAlignment.GetEndPosition() >= Region.LeftPosition ) return WITHIN_REGION;
	    else return BEFORE_REGION;
	}
    }

    // alignment starts on a reference after the left bound
    else {

	// if region has a right boundary
	if ( Region.isRightBoundSpecified() ) {

	    // alignment is on reference between boundaries
	    if ( bAlignment.RefID < Region.RightRefID ) return WITHIN_REGION;

	    // alignment is on reference after right boundary
	    else if ( bAlignment.RefID > Region.RightRefID ) return AFTER_REGION;

	    // alignment is on right bound reference
	    else {
		// check if alignment starts before or at right boundary
		if ( bAlignment.Position <= Region.RightPosition ) return WITHIN_REGION;
		else return AFTER_REGION;
	    }
	}

	// otherwise, alignment is after left bound reference, but there is no right boundary
	else return WITHIN_REGION;
    }
}

// load BAM header data
void BamReaderPrivate::LoadHeaderData(void) {

//    m_header = new BamHeader(&mBGZF);
//    bool headerLoadedOk = m_header->Load();
//    if ( !headerLoadedOk )
//	cerr << "BamReader could not load header" << endl;

    // check to see if proper BAM header
    char buffer[4];
    if (mBGZF.Read(buffer, 4) != 4) {
	fprintf(stderr, "Could not read header type\n");
	exit(1);
    }

    if (strncmp(buffer, "BAM\001", 4)) {
	fprintf(stderr, "wrong header type!\n");
	exit(1);
    }

    // get BAM header text length
    mBGZF.Read(buffer, 4);
    unsigned int headerTextLength = BgzfData::UnpackUnsignedInt(buffer);
    if ( IsBigEndian ) SwapEndian_32(headerTextLength);

    // get BAM header text
    char* headerText = (char*)calloc(headerTextLength + 1, 1);
    mBGZF.Read(headerText, headerTextLength);
    HeaderText = (string)((const char*)headerText);

    // clean up calloc-ed temp variable
    free(headerText);
}

// load existing index data from BAM index file (".bti" OR ".bai"), return success/fail
bool BamReaderPrivate::LoadIndex(const bool lookForIndex, const bool preferStandardIndex) {

    // clear out any existing index data
    ClearIndex();

    // if no index filename provided, so we need to look for available index files
    if ( IndexFilename.empty() ) {

	// attempt to load BamIndex based on current Filename provided & preferStandardIndex flag
	const BamIndex::PreferredIndexType type = (preferStandardIndex ? BamIndex::STANDARD : BamIndex::BAMTOOLS);
	Index = BamIndex::FromBamFilename(Filename, &mBGZF, Parent, type);

	// if null, return failure
	if ( Index == 0 ) return false;

	// generate proper IndexFilename based on type of index created
	IndexFilename = Filename + Index->Extension();
    }

    else {

	// attempt to load BamIndex based on IndexFilename provided by client
	Index = BamIndex::FromIndexFilename(IndexFilename, &mBGZF, Parent);

	// if null, return failure
	if ( Index == 0 ) return false;
    }

    // set cache mode for BamIndex
    Index->SetCacheMode(IndexCacheMode);

    // loading the index data from file
    HasIndex = Index->Load(IndexFilename);

    // mark empty references
    MarkReferences();

    // return index status
    return HasIndex;
}

// populates BamAlignment with alignment data under file pointer, returns success/fail
bool BamReaderPrivate::LoadNextAlignment(BamAlignment& bAlignment) {

    // read in the 'block length' value, make sure it's not zero
    char buffer[4];
    mBGZF.Read(buffer, 4);
    bAlignment.SupportData.BlockLength = BgzfData::UnpackUnsignedInt(buffer);
    if ( IsBigEndian ) { SwapEndian_32(bAlignment.SupportData.BlockLength); }
    if ( bAlignment.SupportData.BlockLength == 0 ) return false;

    // read in core alignment data, make sure the right size of data was read
    char x[BAM_CORE_SIZE];
    if ( mBGZF.Read(x, BAM_CORE_SIZE) != BAM_CORE_SIZE ) return false;

    if ( IsBigEndian ) {
	for ( int i = 0; i < BAM_CORE_SIZE; i+=sizeof(uint32_t) )
	    SwapEndian_32p(&x[i]);
    }

    // set BamAlignment 'core' and 'support' data
    bAlignment.RefID    = BgzfData::UnpackSignedInt(&x[0]);
    bAlignment.Position = BgzfData::UnpackSignedInt(&x[4]);

    unsigned int tempValue = BgzfData::UnpackUnsignedInt(&x[8]);
    bAlignment.Bin        = tempValue >> 16;
    bAlignment.MapQuality = tempValue >> 8 & 0xff;
    bAlignment.SupportData.QueryNameLength = tempValue & 0xff;

    tempValue = BgzfData::UnpackUnsignedInt(&x[12]);
    bAlignment.AlignmentFlag = tempValue >> 16;
    bAlignment.SupportData.NumCigarOperations = tempValue & 0xffff;

    bAlignment.SupportData.QuerySequenceLength = BgzfData::UnpackUnsignedInt(&x[16]);
    bAlignment.MateRefID    = BgzfData::UnpackSignedInt(&x[20]);
    bAlignment.MatePosition = BgzfData::UnpackSignedInt(&x[24]);
    bAlignment.InsertSize   = BgzfData::UnpackSignedInt(&x[28]);

    // set BamAlignment length
    bAlignment.Length = bAlignment.SupportData.QuerySequenceLength;

    // read in character data - make sure proper data size was read
    bool readCharDataOK = false;
    const unsigned int dataLength = bAlignment.SupportData.BlockLength - BAM_CORE_SIZE;
    char* allCharData = (char*)calloc(sizeof(char), dataLength);

    if ( mBGZF.Read(allCharData, dataLength) == (signed int)dataLength) {

	// store 'allCharData' in supportData structure
	bAlignment.SupportData.AllCharData.assign((const char*)allCharData, dataLength);

	// set success flag
	readCharDataOK = true;

	// save CIGAR ops
	// need to calculate this here so that  BamAlignment::GetEndPosition() performs correctly,
	// even when GetNextAlignmentCore() is called
	const unsigned int cigarDataOffset = bAlignment.SupportData.QueryNameLength;
	uint32_t* cigarData = (uint32_t*)(allCharData + cigarDataOffset);
	CigarOp op;
	bAlignment.CigarData.clear();
	bAlignment.CigarData.reserve(bAlignment.SupportData.NumCigarOperations);
	for (unsigned int i = 0; i < bAlignment.SupportData.NumCigarOperations; ++i) {

	    // swap if necessary
	    if ( IsBigEndian ) SwapEndian_32(cigarData[i]);

	    // build CigarOp structure
	    op.Length = (cigarData[i] >> BAM_CIGAR_SHIFT);
	    op.Type   = CIGAR_LOOKUP[ (cigarData[i] & BAM_CIGAR_MASK) ];

	    // save CigarOp
	    bAlignment.CigarData.push_back(op);
	}
    }

    free(allCharData);
    return readCharDataOK;
}

// loads reference data from BAM file
void BamReaderPrivate::LoadReferenceData(void) {

    // get number of reference sequences
    char buffer[4];
    mBGZF.Read(buffer, 4);
    unsigned int numberRefSeqs = BgzfData::UnpackUnsignedInt(buffer);
    if ( IsBigEndian ) SwapEndian_32(numberRefSeqs);
    if ( numberRefSeqs == 0 ) return;
    References.reserve((int)numberRefSeqs);

    // iterate over all references in header
    for (unsigned int i = 0; i != numberRefSeqs; ++i) {

	// get length of reference name
	mBGZF.Read(buffer, 4);
	unsigned int refNameLength = BgzfData::UnpackUnsignedInt(buffer);
	if ( IsBigEndian ) SwapEndian_32(refNameLength);
	char* refName = (char*)calloc(refNameLength, 1);

	// get reference name and reference sequence length
	mBGZF.Read(refName, refNameLength);
	mBGZF.Read(buffer, 4);
	int refLength = BgzfData::UnpackSignedInt(buffer);
	if ( IsBigEndian ) SwapEndian_32(refLength);

	// store data for reference
	RefData aReference;
	aReference.RefName   = (string)((const char*)refName);
	aReference.RefLength = refLength;
	References.push_back(aReference);

	// clean up calloc-ed temp variable
	free(refName);
    }
}

// mark references with no alignment data
void BamReaderPrivate::MarkReferences(void) {

    // ensure index is available
    if ( !HasIndex ) return;

    // mark empty references
    for ( int i = 0; i < (int)References.size(); ++i )
	References.at(i).RefHasAlignments = Index->HasAlignments(i);
}

// opens BAM file (and index)
bool BamReaderPrivate::Open(const string& filename, const string& indexFilename, const bool lookForIndex, const bool preferStandardIndex) {

    // store filenames
    Filename = filename;
    IndexFilename = indexFilename;

    // open the BGZF file for reading, return false on failure
    if ( !mBGZF.Open(filename, "rb") ) return false;

    // retrieve header text & reference data
    LoadHeaderData();
    LoadReferenceData();

    // store file offset of first alignment
    AlignmentsBeginOffset = mBGZF.Tell();

    // if no index filename provided
    if ( IndexFilename.empty() ) {

	// client did not specify that index SHOULD be found
	// useful for cases where sequential access is all that is required
	if ( !lookForIndex ) return true;

	// otherwise, look for index file, return success/fail
	return LoadIndex(lookForIndex, preferStandardIndex) ;
    }

    // client supplied an index filename
    // attempt to load index data, return success/fail
    return LoadIndex(lookForIndex, preferStandardIndex);
}

// returns BAM file pointer to beginning of alignment data
bool BamReaderPrivate::Rewind(void) {

    // rewind to first alignment, return false if unable to seek
    if ( !mBGZF.Seek(AlignmentsBeginOffset) ) return false;

    // retrieve first alignment data, return false if unable to read
    BamAlignment al;
    if ( !LoadNextAlignment(al) ) return false;

    // reset default region info using first alignment in file
    Region.clear();
    HasAlignmentsInRegion = true;

    // rewind back to beginning of first alignment
    // return success/fail of seek
    return mBGZF.Seek(AlignmentsBeginOffset);
}

// change the index caching behavior
void BamReaderPrivate::SetIndexCacheMode(const BamIndex::BamIndexCacheMode mode) {
    IndexCacheMode = mode;
    if ( Index == 0 ) return;
    Index->SetCacheMode(mode);
}

// asks Index to attempt a Jump() to specified region
// returns success/failure
bool BamReaderPrivate::SetRegion(const BamRegion& region) {

    // clear out any prior BamReader region data
    //
    // N.B. - this is cleared so that BamIndex now has free reign to call
    // GetNextAlignmentCore() and do overlap checking without worrying about BamReader
    // performing any overlap checking of its own and moving on to the next read... Calls
    // to GetNextAlignmentCore() with no Region set, simply return the next alignment.
    // This ensures that the Index is able to do just that. (All without exposing
    // LoadNextAlignment() to the public API, and potentially confusing clients with the nomenclature)
    Region.clear();

    // check for existing index
    if ( !HasIndex ) return false;

    // adjust region if necessary to reflect where data actually begins
    BamRegion adjustedRegion(region);
    AdjustRegion(adjustedRegion);

    // if no data present, return true
    // not an error, but BamReader knows that no data is there for future alignment access
    // (this is useful in a MultiBamReader setting where some BAM files may lack data in regions
    // that other BAMs have data)
    if ( !HasAlignmentsInRegion ) {
	Region = adjustedRegion;
	return true;
    }

    // attempt jump to user-specified region return false if jump could not be performed at all
    // (invalid index, unknown reference, etc)
    //
    // Index::Jump() is allowed to modify the HasAlignmentsInRegion flag
    //  * This covers case where a region is requested that lies beyond the last alignment on a reference
    //    If this occurs, any subsequent calls to GetNexAlignment[Core] simply return false
    //    BamMultiReader is then able to successfully pull alignments from a region from multiple files
    //    even if one or more have no data.
    if ( !Index->Jump(adjustedRegion, &HasAlignmentsInRegion) ) return false;

    // save region and return success
    Region = adjustedRegion;
    return true;
}
