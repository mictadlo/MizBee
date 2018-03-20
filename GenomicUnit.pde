
class GenomicUnit   {
  long _length, _start, _stop; // length in bp
  float _rel_start, _rel_stop; // [0,1] length relative to topmost unit total length
  float _container_rel_start, _container_rel_stop;
  boolean _orientation; // T: in increasing coord direction
  int _ref_number; // chromosome reference number [0,n-1]
  float _start_degree, _middle_degree, _stop_degree;
  float[] _render1, _render2, _render3;
  float _or_start_degree, _or_middle_degree, _or_stop_degree;
  float[] _or_render1, _or_render2;
  color _color;
  String _tag;

  ArrayList _subunit;
  int _num_subunits;

  float _pair_value;
  boolean _pair_orientation; // T: orientation is a match, F: inversion

    GenomicUnit( long start, long stop, long len_container, int ref_num ) { 
    _length = stop - start;
    _start = start;
    _stop = stop;
    
    Long l_val1, l_val2;
    l_val1 = start-1;
    l_val2 = len_container-1;
    _rel_start = l_val1.floatValue()/l_val2.floatValue();
    l_val1 = stop-1;
    _rel_stop = l_val1.floatValue()/l_val2.floatValue();
    
    _ref_number = ref_num;

    _orientation = true;

    _subunit = new ArrayList();
    _num_subunits = 0;

    _start_degree = _middle_degree = _stop_degree = 0.0;
    _or_start_degree = _or_middle_degree = _or_stop_degree = 0.0;

    _tag = "";
  }

  void setContainerRelStartAndStop( long container_len, long container_start ) {
    Long l_val1, l_val2;
    l_val1 = _start - container_start - 1;
    l_val2 = container_len-1;
    _container_rel_start = l_val1.floatValue()/l_val2.floatValue();
    l_val1 = _stop - container_start - 1;
    _container_rel_stop = l_val1.floatValue()/l_val2.floatValue();
  }

  float containerRelStart() {
    return _container_rel_start;
  }

  float containerRelStop() {
    return _container_rel_stop;
  }

  boolean locationWithin( long location ) {
    if ( (location >= _start) && (location <= _stop) ) {
      return true;
    }
    return false;
  }

  long midPointCoordinate() {
    return (_stop - _length/2);
  }

  void setPolarRenderPositions( float in_r, float out_r ) {
    _middle_degree = _start_degree + (_stop_degree - _start_degree)*0.5;
    setRenderOnePosition( in_r*cos(_middle_degree), in_r*sin(_middle_degree) );
    setRenderTwoPosition( out_r*cos(_middle_degree), out_r*sin(_middle_degree) );
    //setRenderThreePosition( in_r*CURVE_CP_FRACTION*cos(_middle_degree), in_r*CURVE_CP_FRACTION*sin(_middle_degree) );
  }

  void setStartAndStopDegree( float start_deg, float stop_deg ) {
    _start_degree = start_deg;
    _stop_degree = stop_deg;
  }

  void setMiddleDegree( float mid_deg ) {
    _middle_degree = mid_deg;
  }    

  float startDegree() {
    return _start_degree;
  }

  float stopDegree() {
    return _stop_degree;
  }

  float middleDegree() {
    return _middle_degree;
  }
  
  float halfDegreeLength() {
    return _middle_degree - _start_degree;
  }
  
  void setOuterRingPolarRenderPositions( float in_r, float out_r ) {
    _or_middle_degree = _or_start_degree + (_or_stop_degree - _or_start_degree)*0.5;
    setOuterRingRenderOnePosition( in_r*cos(_or_middle_degree), in_r*sin(_or_middle_degree) );
    setOuterRingRenderTwoPosition( out_r*cos(_or_middle_degree), out_r*sin(_or_middle_degree) );
  }
  
  void setOuterRingStartAndStopDegree( float start_deg, float stop_deg ) {
    _or_start_degree = start_deg;
    _or_stop_degree = stop_deg;
  }

  void setOuterRingMiddleDegree( float mid_deg ) {
    _or_middle_degree = mid_deg;
  }    

  float outerRingStartDegree() {
    return _or_start_degree;
  }

  float outerRingStopDegree() {
    return _or_stop_degree;
  }

  float outerRingMiddleDegree() {
    return _or_middle_degree;
  }
  
  float outerRingHalfDegreeLength() {
    return _or_middle_degree - _or_start_degree;
  }

  void setRenderOnePosition( float x, float y ) {
    _render1 = new float[2];
    _render1[0] = x;
    _render1[1] = y;
  }

  float[] renderOnePosition() {
    return _render1;
  }

  void setRenderTwoPosition( float x, float y ) {
    _render2 = new float[2];
    _render2[0] = x;
    _render2[1] = y;
  }

  float[] renderTwoPosition() {
    return _render2;
  }
  
  void setRenderThreePosition( float x, float y ) {
    _render3 = new float[2];
    _render3[0] = x;
    _render3[1] = y;
  }

  float[] renderThreePosition() {
    return _render3;
  }
  
  void setOuterRingRenderOnePosition( float x, float y ) {
    _or_render1 = new float[2];
    _or_render1[0] = x;
    _or_render1[1] = y;
  }

  float[] outerRingRenderOnePosition() {
    return _or_render1;
  }

  void setOuterRingRenderTwoPosition( float x, float y ) {
    _or_render2 = new float[2];
    _or_render2[0] = x;
    _or_render2[1] = y;
  }

  float[] outerRingRenderTwoPosition() {
    return _or_render2;
  }

  void setColor( color c ) {
    _color = c;
  }

  color getColor() {
    return _color;
  }

  void setOrientation( boolean o ) {
    _orientation = o;
  }

  int refNumber() {
    return _ref_number;
  }

  long start() {
    return _start;
  }

  long stop() {
    return _stop;
  }

  float relStart() {
    return _rel_start;
  }

  float relStop() {
    return _rel_stop;
  }

  long len() {
    return _length;
  }

  boolean orientation() {
    return _orientation;
  }

  boolean empty() {
    if ( _num_subunits == 0 ) {
      return true;
    }
    else {
      return false;
    }
  }

  void trim() {
    _subunit.trimToSize();
    for ( int i = 0; i < _num_subunits; i++ ) {
      ((GenomicUnit)_subunit.get(i)).trim();
    }
  }

  int numSubunits() {
    return _num_subunits;
  }

  void setPairValue( float pv ) {
    _pair_value = pv;
  }

  float pairValue() {
    return _pair_value;
  }

  void setPairOrientation( boolean po ) {
    _pair_orientation = po;
  }

  boolean pairOrientation() {
    return _pair_orientation;
  }

  void setTag( String t ) {
    _tag = t;
  }

  String tag() {
    return _tag;
  }
  
    
  

}

/*****************************************************/

class Chromosome extends GenomicUnit {  
  ArrayList[] _tracks;
  float _number_x, _number_y, _or_number_x, _or_number_y;

  Chromosome( long start, long stop, long len_container, int ref_num ) {
    super( start, stop, len_container, ref_num );

  }

  void addBlock( Block r ) {
    _subunit.add( r );
    ++_num_subunits;

    r.setContainerRelStartAndStop( _length, _start );
  }

  void removeBlock( int index ) {
    _subunit.remove( index );
    --_num_subunits;
  }

  Block getBlock( int index ) {
    return (Block)_subunit.get(index);
  }
  
  void setNumTracks( int num ) {
    _tracks = new ArrayList[num];
    for ( int i = 0; i < num; i++ ) _tracks[i] = new ArrayList();
  }
  
  void addAnnotation( int track, Annotation a ) {
    _tracks[track].add( a );
  }
  
  Annotation getAnnotation( int track, int index ) {
    return (Annotation)_tracks[track].get( index );
  }
  
  int numAnnotations( int track ) {
    return _tracks[track].size();
  }

  void setNumberXAndY( float x, float y ) {
    _number_x = x;
    _number_y = y;
  }

  float numberX() {
    return _number_x;
  }

  float numberY() {
    return _number_y;
  }
  
  void setOuterRingNumberXAndY( float x, float y ) {
    _or_number_x = x;
    _or_number_y = y;
  }

  float outerRingNumberX() {
    return _or_number_x;
  }

  float outerRingNumberY() {
    return _or_number_y;
  }

  float convertDegreeToRelativePosition( float deg ) {  
    return ((deg - _start_degree) / (_stop_degree - _start_degree)); 
  }

  float convertRelativePositionToDegree( float rel_pos ) {
    return ( _start_degree + rel_pos*(_stop_degree - _start_degree) );
  }

}

/*****************************************************/


class Block extends GenomicUnit implements Comparable {
  Block _pair;
  float _min_gene_block_ratio;
  ArrayList _cp;

  Block( long start, long stop, long len_container, int ref_num ) {
    super( start, stop, len_container, ref_num );
    _pair = null;
    
    _cp = new ArrayList();
  }
  
  void clearCP() {
    _cp.clear();
  }
  
  void addCP( float[] pos ) {
    float[] p = new float[2];
    p[0] = pos[0];
    p[1] = pos[1];
    _cp.add( p );
  }
  
  float[] getCP( int index ) {
    return (float[])_cp.get(index);
  }
  
  int numCP() {
    return _cp.size();
  }

  void addGene( Gene g ) {
    _subunit.add( g );
    ++_num_subunits;

    g.setContainerRelStartAndStop( _length, _start );
  }

  void removeGene( int index ) {
    _subunit.remove( index );
    --_num_subunits;
  }

  Gene getGene( int index ) {
    return (Gene)_subunit.get(index);
  }

  void setPair( Block b ) {
    _pair = b;
  }

  Block pair() {
    return _pair;
  }
  
  float minGeneBlockRatio() {
    return _min_gene_block_ratio;
  }
  
  int compareTo(Object another_c)  {
    Long d1 = new Long(_start);
    Long d2 = new Long(((Block) another_c)._start);
    return d1.compareTo(d2);    
  }

}

/*****************************************************/

class Gene extends GenomicUnit {
  boolean _paired;
  Gene _pair;
  int _id;

  Gene( long start, long stop, long len_container, int ref_num ) {
    super( start, stop, len_container, ref_num );
    _pair = null;
  }

  void setPair( Gene g ) {
    _pair = g;
  }

  Gene pair() {
    return _pair;
  }  

  void setID( int id ) {
    _id = id;
  }

  int id() {
    return _id;
  }

}

/*****************************************************/

class Annotation extends GenomicUnit {

  Annotation( long start, long stop, long len_container, int ref_num ) {
    super( start, stop, len_container, ref_num );
  } 
}



