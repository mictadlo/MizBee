
float EXTRA_SPACE_FRACTION = 0.1;
float OR_IN_RADIUS_FRACTION = 0.9;
float IN_RADIUS_FRACTION = 0.7;
float TEXT_OFFSET = 20.0;
float TEXT_DIAMETER = 1.75*TEXT_OFFSET;
float DEG_STEPSIZE = 1.0*DEG_TO_RAD;

color _chrom_color = color( CHROM_COLOR );
float CURVE_CP_FRACTION = 0.55;

int _selected_chromosome = 0;



/*class Combo implements Comparable {
  int _i, _j;
  color _c;
  int _num_elements;
  
  Combo( int i, int j, int num_elements ) {
    _i = i;
    _j = j;
    _num_elements = num_elements;
  }
  
  public int compareTo(Object another_c)  {
    Integer d1 = new Integer(_num_elements);
    Integer d2 = new Integer(((Combo)another_c)._num_elements);
    return d2.compareTo(d1);    
  }
}*/

class GenomeView
{  
  Genome _genome;
  Window _window;
  int _width;
  float _out_radius, _in_radius, _deg_inbetween;
  float _out_diameter, _in_diameter;
  float _or_out_radius, _or_in_radius, _or_out_diameter, _or_in_diameter, _or_deg_inbetween;
  float[] _center;
  
  float _alpha;
  int _alpha_slider_w, _alpha_slider_r;
  int _alpha_slider_x, _alpha_slider_y;
  
  GenomeView( Genome g, int w )
  {
    _genome = g;    
    _center = new float[2];
    setDimensions( w );
    
    _alpha = 0.4;
    _alpha_slider_r = 10;
    
    setChromosomeBlockDegrees();
    setOuterRingChromosomeBlockDegrees();
    buildColorTable();
    
    // create the window
    _window = new Window( _genome.getChromosome(0).startDegree(), _genome.getChromosome(0).stopDegree(), _in_radius, (_out_radius-_in_radius) );
  }
  
  void setDimensions( int w ) {
    _width = w;
    
    _or_out_radius = _width/2.0 - (TEXT_OFFSET+TEXT_DIAMETER/2.0);
    _or_in_radius = OR_IN_RADIUS_FRACTION * _or_out_radius;
    
    _out_radius = _or_in_radius - 2.0*(TEXT_OFFSET+TEXT_DIAMETER/2.0);
    _in_radius = IN_RADIUS_FRACTION * _out_radius;
    
    _or_out_diameter = 2.0*_or_out_radius;
    _or_in_diameter = 2.0*_or_in_radius;
    _out_diameter = 2.0*_out_radius;
    _in_diameter = 2.0*_in_radius;
    
    _center[0] = _center[1] = (float)_width*0.5;
    
    _alpha_slider_x = 20;
    _alpha_slider_y = w-20;
    _alpha_slider_x -= _center[0];
    _alpha_slider_y -= _center[1];
    _alpha_slider_w = w/5;
  }
  
  void render() {
    textFont( sans_12 );
    textAlign( LEFT, TOP );
    fill( 0 );
    text( "source: " + _genome._genome_tag[0], 0, 0 );
    text( "destination: " + _genome._genome_tag[1], 0, 13 );
    
    
    pushMatrix();
    translate( _center[0], _center[1] );
    
    //smooth();
    
    // 
    // render the outer ring
    //
    
    // fill in the color of the unselected chromosome ring
    noStroke();    
    //fill( _fill_color );
    fill( 255 );
    ellipse( 0.0, 0.0, _or_out_diameter, _or_out_diameter );    
    //fill( _background );
    //ellipse( 0.0, 0.0, _or_in_diameter, _or_in_diameter );
    
    //renderOuterRingSelectedChromosome();    
    renderOuterRingSyntenicBlocks();
    
    // outline the chromosome ring
    noFill();
    stroke( _chrom_color );
    strokeWeight( 1.0 );
    ellipse( 0.0, 0.0, _or_out_diameter, _or_out_diameter );
    fill( _background );
    ellipse( 0.0, 0.0, _or_in_diameter, _or_in_diameter );
    
    renderOuterRingSpaces();   
    renderOuterRingNumbers();
    
    //
    // render the inner ring
    //
    
    // fill in the color of the unselected chromosome ring
    noStroke();    
    fill( _fill_color ); 
    ellipse( 0.0, 0.0, _out_diameter, _out_diameter );  
   
    fill( _background );
    ellipse( 0.0, 0.0, _in_diameter, _in_diameter );
    
    // render the selected chromosome arc
    _window.render();
    
    renderSyntenicBlocks();
    
    // outline the chromosome ring
    noFill();
    stroke( _chrom_color );
    strokeWeight( 1.0 );
    ellipse( 0.0, 0.0, _out_diameter, _out_diameter );
    ellipse( 0.0, 0.0, _in_diameter, _in_diameter );
    
    renderSpaces();   
    renderNumbers();
    
    //renderAlphaSlider();
   
    renderOutlineOfSelected();
    popMatrix();
    
    
  }
  
  void renderOutlineOfSelected() {
    float deg_start = (_genome.getChromosome(_selected_chromosome).outerRingStartDegree());
    float deg_end = (_genome.getChromosome(_selected_chromosome).outerRingStopDegree()); 
    deg_start += TWO_PI;
    deg_end += TWO_PI; 
    
    strokeWeight( 2 );
    noFill();
    stroke( 0 );
    
    arc( 0, 0, _or_in_radius*2.0, _or_in_radius*2.0, deg_start, deg_end );
    arc( 0, 0, _or_out_radius*2.0, _or_out_radius*2.0, deg_start, deg_end );
    
    float x1 = _or_in_radius*cos( deg_start );
    float x2 = _or_out_radius*cos( deg_start );
    float y1 = _or_in_radius*sin( deg_start );
    float y2 = _or_out_radius*sin( deg_start );
    line( x1, y1, x2, y2 );
    
    x1 = _or_in_radius*cos( deg_end );
    x2 = _or_out_radius*cos( deg_end );
    y1 = _or_in_radius*sin( deg_end );
    y2 = _or_out_radius*sin( deg_end );
    line( x1, y1, x2, y2 );
    
    deg_start = (_genome.getChromosome(_selected_chromosome).startDegree());
    deg_end = (_genome.getChromosome(_selected_chromosome).stopDegree()); 
    deg_start += TWO_PI;
    deg_end += TWO_PI; 
    
    arc( 0, 0, _in_radius*2.0, _in_radius*2.0, deg_start, deg_end );
    arc( 0, 0, _out_radius*2.0, _out_radius*2.0, deg_start, deg_end );
    
    x1 = _in_radius*cos( deg_start );
    x2 = _out_radius*cos( deg_start );
    y1 = _in_radius*sin( deg_start );
    y2 = _out_radius*sin( deg_start );
    line( x1, y1, x2, y2 );
    
    x1 = _in_radius*cos( deg_end );
    x2 = _out_radius*cos( deg_end );
    y1 = _in_radius*sin( deg_end );
    y2 = _out_radius*sin( deg_end );
    line( x1, y1, x2, y2 );
    
  }
  
  boolean overAlphaSlider( int x, int y ) {
    if ( (x >= _alpha_slider_x) && (x <= _alpha_slider_x+_alpha_slider_w) && 
    (y >= _alpha_slider_y-_alpha_slider_r) && ( y <= _alpha_slider_y+_alpha_slider_r) ) {
      _alpha = float(x - _alpha_slider_x) / (float)_alpha_slider_w;  
      return true;
    }    
    return false;
  }
  
  void renderAlphaSlider() {
    strokeWeight( 1 );
    stroke( 0 );
    
    line( _alpha_slider_x, _alpha_slider_y-_alpha_slider_r, _alpha_slider_x, _alpha_slider_y+_alpha_slider_r );
    line( _alpha_slider_x+_alpha_slider_w, _alpha_slider_y-_alpha_slider_r, _alpha_slider_x+_alpha_slider_w, _alpha_slider_y+_alpha_slider_r );
    line( _alpha_slider_x, _alpha_slider_y, _alpha_slider_x+_alpha_slider_w, _alpha_slider_y );
    
    stroke( 0 );
    fill( _button_color );
    ellipse( _alpha_slider_x+int(_alpha*_alpha_slider_w), _alpha_slider_y, _alpha_slider_r, _alpha_slider_r );
    
    fill( 0 );
    textAlign( LEFT, BOTTOM );
    textFont( sans_15 );
    text( "saturation", _alpha_slider_x, _alpha_slider_y - (_alpha_slider_r+2) );
    text( "line", _alpha_slider_x, _alpha_slider_y - (_alpha_slider_r+13) );
    
    textAlign( RIGHT, CENTER );
    text( "-", _alpha_slider_x-3, _alpha_slider_y );
    textAlign( LEFT, CENTER ); 
    text( "+", _alpha_slider_x+_alpha_slider_w+3, _alpha_slider_y );
  }
  
  void renderOuterRingSelectedChromosome() {
    //smooth();
    noStroke();
    fill( 255 );  
 
    float deg_start = _genome.getChromosome(_selected_chromosome).outerRingStartDegree();
    float deg_end = _genome.getChromosome(_selected_chromosome).outerRingStopDegree();  

    if ( deg_start < deg_end ) {
      beginShape( QUAD_STRIP );
      {
        for ( float a = deg_start; a < deg_end; a += DEG_STEPSIZE )
        {
          vertex( _or_in_radius*cos(a), _or_in_radius*sin(a) );
          vertex( _or_out_radius*cos(a), _or_out_radius*sin(a) );
        }

        vertex( _or_in_radius*cos(deg_end), _or_in_radius*sin(deg_end) );
        vertex( _or_out_radius*cos(deg_end), _or_out_radius*sin(deg_end) );
      }
      endShape();
    }
    else {
      beginShape( QUAD_STRIP );
      {
        for ( float a = deg_start; a < PI; a += DEG_STEPSIZE )
        {
          vertex( _or_in_radius*cos(a), _or_in_radius*sin(a) );
          vertex( _or_out_radius*cos(a), _or_out_radius*sin(a) );
        }
        for ( float a = -PI; a < deg_end; a += DEG_STEPSIZE )
        {
          vertex( _or_in_radius*cos(a), _or_in_radius*sin(a) );
          vertex( _or_out_radius*cos(a), _or_out_radius*sin(a) );
        }

        vertex( _or_in_radius*cos(deg_end), _or_in_radius*sin(deg_end) );
        vertex( _or_out_radius*cos(deg_end), _or_out_radius*sin(deg_end) );
      }
      endShape();
    }
  }
  
  void renderSpaces() {    
    float offset = 1.0;    
    float start_chrom_deg, next_chrom_deg;
    
    // first render for the selected chromosome
    start_chrom_deg = _genome.getChromosome(_selected_chromosome).stopDegree();
    next_chrom_deg = _genome.getChromosome(_genome.genomeChromStartIndex(1-_selected_genome)).startDegree();
  
    if ( start_chrom_deg > next_chrom_deg ) {
      start_chrom_deg -= TWO_PI; 
    }
    
    // fill in with background color
    noStroke();
    fill( _background );
    beginShape( QUAD_STRIP );   
    for ( float a = start_chrom_deg; a <= next_chrom_deg; a += DEG_STEPSIZE ) {
      vertex( (_in_radius-offset)*cos(a), (_in_radius-offset)*sin(a) );
      vertex( (_out_radius+offset)*cos(a), (_out_radius+offset)*sin(a) );
    }

    vertex( (_in_radius-offset)*cos(next_chrom_deg), (_in_radius-offset)*sin(next_chrom_deg) );
    vertex( (_out_radius+offset)*cos(next_chrom_deg), (_out_radius+offset)*sin(next_chrom_deg) );    
    endShape();
     
    // outline the sides
    stroke( _chrom_color );
    strokeWeight( 1 );
    line( _in_radius*cos(start_chrom_deg), _in_radius*sin(start_chrom_deg),
          _out_radius*cos(start_chrom_deg), _out_radius*sin(start_chrom_deg) );
    line( _in_radius*cos(next_chrom_deg), _in_radius*sin(next_chrom_deg),
          _out_radius*cos(next_chrom_deg), _out_radius*sin(next_chrom_deg) );
            
    for ( int i = _genome.genomeChromStartIndex(1-_selected_genome); i < _genome.genomeChromStartIndex(1-_selected_genome)+_genome.numChromosomes(1-_selected_genome); i++ ) {
      start_chrom_deg = _genome.getChromosome(i).stopDegree();
      if ( i == (_genome.genomeChromStartIndex(1-_selected_genome)+_genome.numChromosomes(1-_selected_genome)-1) ) {
        next_chrom_deg = _genome.getChromosome(_selected_chromosome).startDegree();
      }
      else {
        next_chrom_deg = _genome.getChromosome(i+1).startDegree();
      }
  
      if ( start_chrom_deg > next_chrom_deg ) {
        start_chrom_deg -= TWO_PI; 
      }
    
      // fill in with background color
      noStroke();
      fill( _background );
      beginShape( QUAD_STRIP );   
      for ( float a = start_chrom_deg; a <= next_chrom_deg; a += DEG_STEPSIZE ) {
        vertex( (_in_radius-offset)*cos(a), (_in_radius-offset)*sin(a) );
        vertex( (_out_radius+offset)*cos(a), (_out_radius+offset)*sin(a) );
      }

      vertex( (_in_radius-offset)*cos(next_chrom_deg), (_in_radius-offset)*sin(next_chrom_deg) );
      vertex( (_out_radius+offset)*cos(next_chrom_deg), (_out_radius+offset)*sin(next_chrom_deg) );    
      endShape();
     
      // outline the sides
      stroke( _chrom_color );
      strokeWeight( 1 );
      line( _in_radius*cos(start_chrom_deg), _in_radius*sin(start_chrom_deg),
            _out_radius*cos(start_chrom_deg), _out_radius*sin(start_chrom_deg) );
      line( _in_radius*cos(next_chrom_deg), _in_radius*sin(next_chrom_deg),
            _out_radius*cos(next_chrom_deg), _out_radius*sin(next_chrom_deg) );
    }
  }
  
  void renderOuterRingSpaces() {    
    float offset = 1.0;    
    float start_chrom_deg, next_chrom_deg;
    for ( int i = _genome.genomeChromStartIndex(_selected_genome); i < (_genome.genomeChromStartIndex(_selected_genome)+_genome.numChromosomes(_selected_genome)); i++ ) {
      start_chrom_deg = _genome.getChromosome(i).outerRingStopDegree();
      if ( i == (_genome.genomeChromStartIndex(_selected_genome)+_genome.numChromosomes(_selected_genome)-1) ) {
        next_chrom_deg = _genome.getChromosome(_genome.genomeChromStartIndex(_selected_genome)).outerRingStartDegree();
      }
      else {
        next_chrom_deg = _genome.getChromosome(i+1).outerRingStartDegree();
      }
  
      if ( start_chrom_deg > next_chrom_deg ) {
        start_chrom_deg -= TWO_PI; 
      }
    
      // fill in with background color
      noStroke();
      fill( _background );
      beginShape( QUAD_STRIP );   
      for ( float a = start_chrom_deg; a <= next_chrom_deg; a += DEG_STEPSIZE ) {
        vertex( (_or_in_radius-offset)*cos(a), (_or_in_radius-offset)*sin(a) );
        vertex( (_or_out_radius+offset)*cos(a), (_or_out_radius+offset)*sin(a) );
      }

      vertex( (_or_in_radius-offset)*cos(next_chrom_deg), (_or_in_radius-offset)*sin(next_chrom_deg) );
      vertex( (_or_out_radius+offset)*cos(next_chrom_deg), (_or_out_radius+offset)*sin(next_chrom_deg) );    
      endShape();
     
      // outline the sides
      stroke( _chrom_color );
      strokeWeight( 1 );
      line( _or_in_radius*cos(start_chrom_deg), _or_in_radius*sin(start_chrom_deg),
            _or_out_radius*cos(start_chrom_deg), _or_out_radius*sin(start_chrom_deg) );
      line( _or_in_radius*cos(next_chrom_deg), _or_in_radius*sin(next_chrom_deg),
            _or_out_radius*cos(next_chrom_deg), _or_out_radius*sin(next_chrom_deg) );
    }
  }
  
  void renderNumbers() {
    //smooth();
    noStroke();
    textFont( sans_12 );
    textAlign( CENTER, CENTER );
    
    float num_x, num_y, rot_deg, deg_change = 0.0;
    num_x = _genome.getChromosome(_selected_chromosome).numberX();
    num_y = _genome.getChromosome(_selected_chromosome).numberY();
    rot_deg = _genome.getChromosome(_selected_chromosome).middleDegree();
    if ( (rot_deg < PI) && (rot_deg > 0.0) ) deg_change = -HALF_PI;
    else deg_change = HALF_PI;
    fill( 0 );
    
    float half_word_len = 0.5*textWidth(_genome.getChromosome(_selected_chromosome).tag());
    float half_arc_length_approx = _out_radius * sin( _genome.getChromosome(_selected_chromosome).halfDegreeLength() );
    if ( half_word_len > half_arc_length_approx ) {
      if ( (rot_deg < (HALF_PI)) && (rot_deg > -HALF_PI) ) deg_change = 0.0;
      else deg_change = PI;
    }
    rot_deg += deg_change;
    
    pushMatrix();
    translate( num_x, num_y );
    rotate( rot_deg );
    text( _genome.getChromosome(_selected_chromosome).tag(), 0, 0 );
    popMatrix();
        
    for ( int i = _genome.genomeChromStartIndex(1-_selected_genome); i < _genome.genomeChromStartIndex(1-_selected_genome)+_genome.numChromosomes(1-_selected_genome); i++ ) {
      deg_change = 0.0;
      num_x = _genome.getChromosome(i).numberX();
      num_y = _genome.getChromosome(i).numberY();
      rot_deg = _genome.getChromosome(i).middleDegree();
      if ( (rot_deg < PI) && (rot_deg > 0.0) ) deg_change = -HALF_PI;
      else deg_change = HALF_PI;
      
      fill( _chrom_color );
      
      half_word_len = 0.5*textWidth(_genome.getChromosome(i).tag());
      half_arc_length_approx = _out_radius * sin( _genome.getChromosome(i).halfDegreeLength() );
      if ( half_word_len > half_arc_length_approx ) {
        if ( (rot_deg < (HALF_PI)) && (rot_deg > -HALF_PI) ) deg_change = 0.0;
        else deg_change = PI;
      }
      rot_deg += deg_change;
    
      
      pushMatrix();
      translate( num_x, num_y );
      rotate( rot_deg );
      text( _genome.getChromosome(i).tag(), 0, 0 );
      popMatrix();
    }  
  }
  
  void renderOuterRingNumbers() {
    //smooth();
    noStroke();
    textFont( sans_12 );
    textAlign( CENTER, CENTER );
    
    float num_x, num_y, rot_deg, deg_change, half_word_len, half_arc_length_approx;
    for ( int i = _genome.genomeChromStartIndex(_selected_genome); i < _genome.genomeChromStartIndex(_selected_genome)+_genome.numChromosomes(_selected_genome); i++ ) {
      deg_change = 0.0;
      num_x = _genome.getChromosome(i).outerRingNumberX();
      num_y = _genome.getChromosome(i).outerRingNumberY();
      
      rot_deg = _genome.getChromosome(i).outerRingMiddleDegree();
      if ( (rot_deg < PI) && (rot_deg > 0.0) ) deg_change = -HALF_PI;
      else deg_change = HALF_PI;
      
      half_word_len = 0.5*textWidth(_genome.getChromosome(i).tag());
      half_arc_length_approx = _or_out_radius * sin( _genome.getChromosome(i).outerRingHalfDegreeLength() );
      if ( half_word_len > half_arc_length_approx ) {
        if ( (rot_deg < (HALF_PI)) && (rot_deg > -HALF_PI) ) deg_change = 0.0;
        else deg_change = PI;
      }
      rot_deg += deg_change;
      
      pushMatrix();
      translate( num_x, num_y );
      rotate( rot_deg );
      
      if ( i == _selected_chromosome ) {
        fill( 255 );
        ellipse( 0, 0, TEXT_DIAMETER, TEXT_DIAMETER );
        fill( 0 );
      }
      else {
        fill( _chrom_color );
      }
      text( _genome.getChromosome(i).tag(), 0, 0 );
      
      popMatrix();
    }
  }
  
  void renderSyntenicBlocks() {
    float shutter_start = _window.startDegree();
    float shutter_end = _window.endDegree();
    
    float chrom_start, chrom_end, render_start, render_end;
    // baiscally always doing this part of the if statement....
    if ( shutter_end >= shutter_start ) {

        // most of this stuff is for dealing with the "window"
        chrom_start = _genome.getChromosome(_selected_chromosome).startDegree();
        chrom_end = _genome.getChromosome(_selected_chromosome).stopDegree();
        
        if ( (chrom_end > shutter_start) || (chrom_start < shutter_end) ) {
          render_start = max( chrom_start, shutter_start );
          render_end = min( chrom_end, shutter_end );
          
          // you will call this!!!!!
          renderChromosomeBlocks( _selected_chromosome, render_start, render_end );
        }

    }
    else {
        chrom_start = _genome.getChromosome(_selected_chromosome).startDegree();
        chrom_end = _genome.getChromosome(_selected_chromosome).stopDegree();
        
        if ( chrom_start <= shutter_end ) {
          render_start = chrom_start;
          render_end = min( chrom_end, shutter_end );
          
          renderChromosomeBlocks( _selected_chromosome, render_start, render_end );
        }
        if ( chrom_end > shutter_start ) {
          render_start = max( chrom_start, shutter_start );
          render_end = chrom_end;
          
          renderChromosomeBlocks( _selected_chromosome, render_start, render_end );
        }
    }
  }
  
  void renderOuterRingSyntenicBlocks() {    
    float[] s_pos1, e_pos1;
    s_pos1 = new float[2];
    e_pos1 = new float[2];
    //smooth();
    strokeWeight( 1 );
    for ( int i = _genome.genomeChromStartIndex(_selected_genome); i < _genome.genomeChromStartIndex(_selected_genome)+_genome.numChromosomes(_selected_genome); i++ ) {
      //stroke( _fill_color );
      for ( int j = 0; j < _genome.getChromosome(i).numSubunits(); j++ ) {
        if ( i == _selected_chromosome ) {
          stroke( _genome.getChromosome(i).getBlock(j).getColor(), _alpha*255.0 );
        }
        else {
          stroke( _genome.getChromosome(i).getBlock(j).getColor(), _alpha*200.0 );
        }
      
        // render this syntenic block
        s_pos1 = _genome.getChromosome(i).getBlock(j).outerRingRenderOnePosition();
        e_pos1 = _genome.getChromosome(i).getBlock(j).outerRingRenderTwoPosition();
                            
        line( s_pos1[0], s_pos1[1], e_pos1[0], e_pos1[1] );
      }
    }
  }
  
  boolean blockSelected( float x, float y ) {

    // convert position into polar coordinates
    float r = sqrt( x*x + y*y );
        
    // check if the r is inbetween the radius values
    if ( (r < _in_radius) || (r > _out_radius) ) {
      //_block_selected = false;
      //_selected_block = 0;
      return false;
    }
    
    // determine if the degree is within the window
    float shutter_start = _window.startDegree();
    float shutter_end = _window.endDegree();
    float theta = (float)Math.atan2( (double)y,(double)(x+1.0e-6) );
    
    if ( shutter_end > shutter_start ) {
      if ( !( (theta >= shutter_start) && (theta <= shutter_end) ) ) {
        //_block_selected = false;
        //_selected_block = 0;
        return false;
      }
    }
    else {
      if ( (theta > shutter_end) && (theta < shutter_start) ) {
        //_block_selected = false;
        //_selected_block = 0;
        return false;
      }
    }
    
    // determine which chromosome we have selected
    Chromosome c;
    float block_start, block_end;
    for ( int i = 0; i < _genome.numChromosomes(0)+_genome.numChromosomes(1); i++ ) {
      c = _genome.getChromosome(i);
           
      if ( (theta >= c.startDegree()) && (theta <= c.stopDegree()) ) {
        // go through the regions and determine which one, if any, we have selected
        for ( int j = 0; j < c.numSubunits(); j++ ) {
          block_start = c.getBlock(j).startDegree();
          block_end = c.getBlock(j).stopDegree();
          
          if ( (theta >= block_start) && (theta <= block_end) ) {
            _selected_chromosome = i;
            _selected_block = j;
            _block_selected = true;
            return true;
          }     
        }
        
        // if no regions were selected, just return false
        //_block_selected = false;
        //_selected_block = 0;
        return false;
      }
    }
    
    //_block_selected = false;
    //_selected_block = 0;
    return false;
  }
  
  void renderChromosomeBlocks( int index, float start_deg, float stop_deg ) {
    float block_start, block_end, block_mid;
    float[] s_pos1, e_pos1, s_pos2, e_pos2, cp_pos1, cp_pos2, cp_pos3, cp_pos4;
    s_pos1 = new float[2];
    e_pos1 = new float[2];
    s_pos2 = new float[2];
    e_pos2 = new float[2];
    cp_pos1 = new float[2];
    cp_pos2 = new float[2];
    cp_pos3 = new float[2];
    cp_pos4 = new float[2];
    for ( int i = 0; i < _genome.getChromosome(index).numSubunits(); i++ ) {
      if ( _block_selected && (index == _selected_chromosome) && (i == _selected_block) ) {
        strokeWeight( 2.0 );
        stroke( color(0) );
      }
      else {
        strokeWeight( 1.0 );
        stroke( _genome.getChromosome(index).getBlock(i).getColor(), _alpha*255.0 );
      }
      
      // checking if block is within "window"
        block_start = _genome.getChromosome(index).getBlock(i).startDegree();
        block_end = _genome.getChromosome(index).getBlock(i).stopDegree();
        
        
        if ( (block_start < stop_deg) && (block_end > start_deg)) {       
          
          // render this syntenic block
          s_pos1 = _genome.getChromosome(index).getBlock(i).renderOnePosition();
          e_pos1 = _genome.getChromosome(index).getBlock(i).renderTwoPosition();
                            
          line( s_pos1[0], s_pos1[1], e_pos1[0], e_pos1[1] );
          
          s_pos2 = _genome.getChromosome(index).getBlock(i).pair().renderOnePosition();
          e_pos2 = _genome.getChromosome(index).getBlock(i).pair().renderTwoPosition();
          
          line( s_pos2[0], s_pos2[1], e_pos2[0], e_pos2[1] );
          
          /*cp_pos1 = _genome.getChromosome(index).getBlock(i).renderThreePosition();
          cp_pos2 = _genome.getChromosome(index).getBlock(i).pair().renderThreePosition();
          
          noFill();  
          bezier( s_pos1[0], s_pos1[1], cp_pos1[0], cp_pos1[1], cp_pos2[0], cp_pos2[1], s_pos2[0], s_pos2[1] );*/
          
          noFill();  
          //bsplineDetail( 1 );
          bsplineTightness( 0.95 );
          
          beginBspline();
          bsplineVertex( s_pos1[0], s_pos1[1] );
          for ( int cp = 0; cp < _genome.getChromosome(index).getBlock(i).numCP(); cp++ ) {
            cp_pos1 = _genome.getChromosome(index).getBlock(i).getCP(cp);
            bsplineVertex( cp_pos1[0], cp_pos1[1] );
          }
          bsplineVertex( s_pos2[0], s_pos2[1] );
          endBspline();

      }
    }
  }
  
  float windowStartDegree() {
    return _window.startDegree();
  }
  
  float windowEndDegree() {
    return _window.endDegree();
  }
  
  boolean overStartShutter( float x, float y ) {
    return _window.overStartShutter( x, y );
  }
  
  boolean overEndShutter( float x, float y ) {
    return _window.overEndShutter( x, y );
  }
  
  void updateShutter( float x, float y, boolean start_shutter, boolean end_shutter, float start_deg, float stop_deg ) {
    _window.updateShutter( x, y, start_shutter, end_shutter, start_deg, stop_deg );
  }
  
  void setWindowStartAndEnd( float start, float end ) {
    _window.setStartAndEndDegrees( start, end );
  }
  
  boolean overChromosomeNumber( float x, float y ) {
    float d;
    float num_x, num_y;
    for ( int i = _genome.genomeChromStartIndex(_selected_genome); i < _genome.genomeChromStartIndex(_selected_genome)+_genome.numChromosomes(_selected_genome); i++ ) {
      num_x = _genome.getChromosome(i).outerRingNumberX();
      num_y = _genome.getChromosome(i).outerRingNumberY();
      
      d = sqrt( (x-num_x)*(x-num_x) + (y-num_y)*(y-num_y) );
      if ( d <= (TEXT_DIAMETER/2.0) ) {
        _selected_chromosome = i;
        setChromosomeBlockDegrees();
        _window.setStartAndEndDegrees( _genome.getChromosome(i).startDegree(), _genome.getChromosome(i).stopDegree() );
        
        _selected_block = 0;
        //_block_selected = false;
        return true;
      }
    }
    return false;
  }
  
  
  
  float centerX() {
    return _center[0];
  }
  
  float centerY() {
    return _center[1];
  }
  
  void setChromosomeBlockDegrees() {
    
    // want to have whitespace inbetween chromosomes when rendered,
    //  so, need to compute the total length and the 
    //  length inbetween (in degrees)
    
    // total extra whitespace in genome units (bp), as a 
    //  percentage of the total length of genome
    long selected_chromosome_len;
    selected_chromosome_len = _genome.getChromosome(_selected_chromosome).len();
    //selected_chromosome_len = _genome.len(1-_selected_genome);
    
    Long l_val = _genome.len(1-_selected_genome)+selected_chromosome_len;
    float total_extra_space = l_val.floatValue() * EXTRA_SPACE_FRACTION;
    
    // now, compute the total length over which we will be
    //  rendering a circle
    l_val = _genome.len(1-_selected_genome)+selected_chromosome_len;
    float total_length = l_val.floatValue() + total_extra_space;
    
    // find the arclength of the extra white space
    float single_extra_space = (total_extra_space/total_length) /
                               float(_genome.numChromosomes(1-_selected_genome)+1);
    _deg_inbetween = TWO_PI * single_extra_space;
    
    // now compute the relative lengths of each chromosome, ie.
    //  the percentage of the genome length plus extra white space 
    float start_deg = -90.0*DEG_TO_RAD;
    float middle_deg, stop_deg;
    float reg_start, reg_end, block_start, block_end, chrom_deg_diff;
    
    // set the start and stop for the selected chromosome, along with the blocks
    middle_deg = start_deg + (float)selected_chromosome_len/total_length * PI;
    stop_deg = start_deg + (float)selected_chromosome_len/total_length * TWO_PI;
           
    if ( middle_deg > PI ) {
      middle_deg -= TWO_PI;
    }
    if ( stop_deg > PI ) {
      stop_deg -= TWO_PI;
    }
      
    if ( start_deg > stop_deg ) {
      start_deg -= TWO_PI;
    }      
      
    _genome.getChromosome(_selected_chromosome).setStartAndStopDegree( start_deg, stop_deg );
    _genome.getChromosome(_selected_chromosome).setMiddleDegree( middle_deg );
    _genome.getChromosome(_selected_chromosome).setNumberXAndY( (_out_radius+TEXT_OFFSET)*cos(middle_deg), (_out_radius+TEXT_OFFSET)*sin(middle_deg) );
      
    chrom_deg_diff = (float)selected_chromosome_len/total_length * TWO_PI;
    float[] chrom_cp = new float[2];
    
    for ( int j = 0; j < _genome.getChromosome(_selected_chromosome).numSubunits(); j++ ) {        
      reg_start = start_deg + chrom_deg_diff * _genome.getChromosome(_selected_chromosome).getBlock(j).relStart();
      reg_end = start_deg + chrom_deg_diff * _genome.getChromosome(_selected_chromosome).getBlock(j).relStop();
        
      if ( reg_start > PI ) {
        reg_start -= TWO_PI;
      }
      if ( reg_end > PI ) {
        reg_end -= TWO_PI;
      }
        
      _genome.getChromosome(_selected_chromosome).getBlock(j).setStartAndStopDegree( reg_start, reg_end );
      _genome.getChromosome(_selected_chromosome).getBlock(j).setPolarRenderPositions( _in_radius, _out_radius );      
    }
      
    start_deg = stop_deg + _deg_inbetween;
       
    // set the degrees for the unselected genome chromosomes
    for ( int i = _genome.genomeChromStartIndex(1-_selected_genome); i < _genome.genomeChromStartIndex(1-_selected_genome)+_genome.numChromosomes(1-_selected_genome); i++ ) {
      middle_deg = start_deg + (float)_genome.getChromosome(i).len()/total_length * PI;
      stop_deg = start_deg + (float)_genome.getChromosome(i).len()/total_length * TWO_PI;
           
      if ( middle_deg > PI ) {
        middle_deg -= TWO_PI;
      }
      if ( stop_deg > PI ) {
        stop_deg -= TWO_PI;
      }
      
      if ( start_deg > stop_deg ) {
        start_deg -= TWO_PI;
      }      
      
      _genome.getChromosome(i).setStartAndStopDegree( start_deg, stop_deg );
      _genome.getChromosome(i).setMiddleDegree( middle_deg );
      _genome.getChromosome(i).setNumberXAndY( (_out_radius+TEXT_OFFSET)*cos(middle_deg), (_out_radius+TEXT_OFFSET)*sin(middle_deg) );
      
      chrom_deg_diff = (float)_genome.getChromosome(i).len()/total_length * TWO_PI;
      for ( int j = 0; j < _genome.getChromosome(i).numSubunits(); j++ ) {        
        reg_start = start_deg + chrom_deg_diff * _genome.getChromosome(i).getBlock(j).relStart();
        reg_end = start_deg + chrom_deg_diff * _genome.getChromosome(i).getBlock(j).relStop();
        
        if ( reg_start > PI ) {
          reg_start -= TWO_PI;
        }
        if ( reg_end > PI ) {
          reg_end -= TWO_PI;
        }
        
        _genome.getChromosome(i).getBlock(j).setStartAndStopDegree( reg_start, reg_end );
        _genome.getChromosome(i).getBlock(j).setPolarRenderPositions( _in_radius, _out_radius );
      }
      
      start_deg = stop_deg + _deg_inbetween;
    }
    
    // set the degrees for the pairs
    int chrom_index;        
    for ( int j = 0; j < _genome.getChromosome(_selected_chromosome).numSubunits(); j++ ) {  
      chrom_index = _genome.getChromosome(_selected_chromosome).getBlock(j).pair().refNumber();
      start_deg = _genome.getChromosome(chrom_index).startDegree();
      stop_deg = _genome.getChromosome(chrom_index).stopDegree();
      if ( start_deg > stop_deg ) {
        start_deg -= TWO_PI;
      } 
      
      chrom_deg_diff = (float)_genome.getChromosome(chrom_index).len()/total_length * TWO_PI;
  
      reg_start = start_deg + chrom_deg_diff * _genome.getChromosome(_selected_chromosome).getBlock(j).pair().relStart();
      reg_end = start_deg  + chrom_deg_diff * _genome.getChromosome(_selected_chromosome).getBlock(j).pair().relStop();
        
      if ( reg_start > PI ) {
        reg_start -= TWO_PI;
      }
      if ( reg_end > PI ) {
        reg_end -= TWO_PI;
      }
        
      _genome.getChromosome(_selected_chromosome).getBlock(j).pair().setStartAndStopDegree( reg_start, reg_end );
      _genome.getChromosome(_selected_chromosome).getBlock(j).pair().setPolarRenderPositions( _in_radius, _out_radius );      
    }
    
    // set the control points for the blocks
    float[] c1, c2, mid, ochrom_cp;
    c1 = new float[2];
    c2 = new float[2];
    mid = new float[2];
    ochrom_cp = new float[2];
    float r, theta;
    int counter;
    float cp1_frac = 0.75;
    float cp2_frac = 0.95;
    float cp3_frac = 0.9;
    for ( int j = 0; j < _genome.getChromosome(_selected_chromosome).numSubunits(); j++ ) { 
      // go through the list and find out how many neighbors this block has
      counter = 0;
      chrom_index = _genome.getChromosome(_selected_chromosome).getBlock(j).pair().refNumber();
      for ( int n = j+1; n < _genome.getChromosome(_selected_chromosome).numSubunits(); n++ ) {
        if ( _genome.getChromosome(_selected_chromosome).getBlock(n).pair().refNumber() == chrom_index ) ++counter;
        else break;
      }
      
      // get the x and y of the center points of the two chromosomes      
      c1[0] = _in_radius*cos(_genome.getChromosome(_selected_chromosome).middleDegree());
      c1[1] = _in_radius*sin(_genome.getChromosome(_selected_chromosome).middleDegree());
      c2[0] = _in_radius*cos(_genome.getChromosome(chrom_index).middleDegree());
      c2[1] = _in_radius*sin(_genome.getChromosome(chrom_index).middleDegree());
      
      mid[0] = c1[0] + (c2[0] - c1[0])/1.25;
      mid[1] = c1[1] + (c2[1] - c1[1])/1.25;
      
      r = sqrt(mid[0]*mid[0] + mid[1]*mid[1]);
      r *= cp1_frac;
      
      theta = atan2( mid[1], mid[0] );
      
      mid[0] = r*cos(theta);
      mid[1] = r*sin(theta);
        
      theta = _genome.getChromosome(_selected_chromosome).getBlock(j).middleDegree()+_genome.getChromosome(_selected_chromosome).getBlock(j+counter).middleDegree();
      theta *= 0.5;
      chrom_cp[0] = _in_radius*cp2_frac*cos(theta);
      chrom_cp[1] = _in_radius*cp2_frac*sin(theta);
      
      /*theta = _genome.getChromosome(_selected_chromosome).getBlock(j).pair().middleDegree()+_genome.getChromosome(_selected_chromosome).getBlock(j+counter).pair().middleDegree();
      theta *= 0.5;R
      ochrom_cp[0] = _in_radius*0.9*cos(theta);
      ochrom_cp[1] = _in_radius*0.9*sin(theta);*/
      
      _genome.getChromosome(_selected_chromosome).getBlock(j).clearCP();
      
      _genome.getChromosome(_selected_chromosome).getBlock(j).addCP( chrom_cp );      
      _genome.getChromosome(_selected_chromosome).getBlock(j).addCP( mid );     
     
      ochrom_cp[0] = _in_radius*cp3_frac*cos(_genome.getChromosome(_selected_chromosome).getBlock(j).pair().middleDegree());
      ochrom_cp[1] = _in_radius*cp3_frac*sin(_genome.getChromosome(_selected_chromosome).getBlock(j).pair().middleDegree()); 
      _genome.getChromosome(_selected_chromosome).getBlock(j).addCP( ochrom_cp );  
      
      for ( int n = 0; n < counter; n++ ) {
        _genome.getChromosome(_selected_chromosome).getBlock(j+n+1).clearCP();
        
        _genome.getChromosome(_selected_chromosome).getBlock(j+n+1).addCP( chrom_cp );      
        _genome.getChromosome(_selected_chromosome).getBlock(j+n+1).addCP( mid );
      
        ochrom_cp[0] = _in_radius*cp3_frac*cos(_genome.getChromosome(_selected_chromosome).getBlock(j+n+1).pair().middleDegree());
        ochrom_cp[1] = _in_radius*cp3_frac*sin(_genome.getChromosome(_selected_chromosome).getBlock(j+n+1).pair().middleDegree());
        _genome.getChromosome(_selected_chromosome).getBlock(j+n+1).addCP( ochrom_cp );  
      }
      
      j += counter;
    }
  }
  
  void setOuterRingChromosomeBlockDegrees() {  
    // the first genome 
    Long l_val = _genome.len(0);
    float total_extra_space = l_val.floatValue() * EXTRA_SPACE_FRACTION;
    float total_length = l_val.floatValue() + total_extra_space;
    
    float single_extra_space = (total_extra_space/total_length) /
                               float(_genome.numChromosomes(0));
    _or_deg_inbetween = TWO_PI * single_extra_space;
    
    float start_deg = -90.0*DEG_TO_RAD;
    float middle_deg, stop_deg;
    float reg_start, reg_end, block_start, block_end, chrom_deg_diff;
    for ( int i = 0; i < _genome.numChromosomes(0); i++ ) {
      middle_deg = start_deg + (float)_genome.getChromosome(i).len()/total_length * PI;
      stop_deg = start_deg + (float)_genome.getChromosome(i).len()/total_length * TWO_PI;
           
      if ( middle_deg > PI ) {
        middle_deg -= TWO_PI;
      }
      if ( stop_deg > PI ) {
        stop_deg -= TWO_PI;
      }
      
      if ( start_deg > stop_deg ) {
        start_deg -= TWO_PI;
      }      
      
      _genome.getChromosome(i).setOuterRingStartAndStopDegree( start_deg, stop_deg );
      _genome.getChromosome(i).setOuterRingMiddleDegree( middle_deg );
      _genome.getChromosome(i).setOuterRingNumberXAndY( (_or_out_radius+TEXT_OFFSET)*cos(middle_deg), (_or_out_radius+TEXT_OFFSET)*sin(middle_deg) );
      
      chrom_deg_diff = (float)_genome.getChromosome(i).len()/total_length * TWO_PI;
      for ( int j = 0; j < _genome.getChromosome(i).numSubunits(); j++ ) {        
        reg_start = start_deg + chrom_deg_diff * _genome.getChromosome(i).getBlock(j).relStart();
        reg_end = start_deg + chrom_deg_diff * _genome.getChromosome(i).getBlock(j).relStop();
        
        if ( reg_start > PI ) {
          reg_start -= TWO_PI;
        }
        if ( reg_end > PI ) {
          reg_end -= TWO_PI;
        }
        
        _genome.getChromosome(i).getBlock(j).setOuterRingStartAndStopDegree( reg_start, reg_end );
        _genome.getChromosome(i).getBlock(j).setOuterRingPolarRenderPositions( _or_in_radius, _or_out_radius );
      }
      
      start_deg = stop_deg + _or_deg_inbetween;
    }
    
    // the other genome
    l_val = _genome.len(1);
    total_extra_space = l_val.floatValue() * EXTRA_SPACE_FRACTION;
    total_length = l_val.floatValue() + total_extra_space;
    
    single_extra_space = (total_extra_space/total_length) /
                          float(_genome.numChromosomes(1));
    
    start_deg = -90.0*DEG_TO_RAD;
    for ( int i = _genome.genomeChromStartIndex(1); i < _genome.genomeChromStartIndex(1)+_genome.numChromosomes(1); i++ ) {
      middle_deg = start_deg + (float)_genome.getChromosome(i).len()/total_length * PI;
      stop_deg = start_deg + (float)_genome.getChromosome(i).len()/total_length * TWO_PI;
           
      if ( middle_deg > PI ) {
        middle_deg -= TWO_PI;
      }
      if ( stop_deg > PI ) {
        stop_deg -= TWO_PI;
      }
      
      if ( start_deg > stop_deg ) {
        start_deg -= TWO_PI;
      }      
     
      _genome.getChromosome(i).setOuterRingStartAndStopDegree( start_deg, stop_deg );
      _genome.getChromosome(i).setOuterRingMiddleDegree( middle_deg );
      _genome.getChromosome(i).setOuterRingNumberXAndY( (_or_out_radius+TEXT_OFFSET)*cos(middle_deg), (_or_out_radius+TEXT_OFFSET)*sin(middle_deg) );
      
      chrom_deg_diff = (float)_genome.getChromosome(i).len()/total_length * TWO_PI;
      for ( int j = 0; j < _genome.getChromosome(i).numSubunits(); j++ ) {        
        reg_start = start_deg + chrom_deg_diff * _genome.getChromosome(i).getBlock(j).relStart();
        reg_end = start_deg + chrom_deg_diff * _genome.getChromosome(i).getBlock(j).relStop();
        
        if ( reg_start > PI ) {
          reg_start -= TWO_PI;
        }
        if ( reg_end > PI ) {
          reg_end -= TWO_PI;
        }
        
        _genome.getChromosome(i).getBlock(j).setOuterRingStartAndStopDegree( reg_start, reg_end );
        _genome.getChromosome(i).getBlock(j).setOuterRingPolarRenderPositions( _or_in_radius, _or_out_radius );
      }
      
      start_deg = stop_deg + _or_deg_inbetween;
    }
  }
  
  void rightKeyPressed() {
    if ( _selected_chromosome == (_genome.genomeChromStartIndex(_selected_genome)+_genome.numChromosomes(_selected_genome)-1) ) {
      _selected_chromosome = _genome.genomeChromStartIndex(_selected_genome);
    }
    else {
      _selected_chromosome += 1;
    }
    
    setChromosomeBlockDegrees();
    _window.setStartAndEndDegrees( _genome.getChromosome(_selected_chromosome).startDegree(), _genome.getChromosome(_selected_chromosome).stopDegree() );
        
    _selected_block = 0;
    _rolled_over_block = 0;
    
    _gene_selected = _gene_rolled_over = _block_rolled_over = false;
    _selected_gene = _rolled_over_gene = 0;
  }
  
  void leftKeyPressed() {
    if ( _selected_chromosome == _genome.genomeChromStartIndex(_selected_genome) ) {
      _selected_chromosome = _genome.genomeChromStartIndex(_selected_genome)+_genome.numChromosomes(_selected_genome)-1;
    }
    else {
      _selected_chromosome -= 1;
    }
    
    setChromosomeBlockDegrees();
    _window.setStartAndEndDegrees( _genome.getChromosome(_selected_chromosome).startDegree(), _genome.getChromosome(_selected_chromosome).stopDegree() );
        
    _selected_block = 0;
    _rolled_over_block = 0;
    
    _gene_selected = _gene_rolled_over = _block_rolled_over = false;
    _selected_gene = _rolled_over_gene = 0;
  }
  
  void buildColorTable() {
    
    // build color map
    color[] color_map = { color(#e31a1c),
                          color(#377db8),
                          color(#4daf4a),
                          color(#984ea3),
                          color(#ff7f00),
                          color(#ffff33),
                          color(#a65628),
                          color(#f781bf) };
    
    // build color table, ie for each chrom, assign a color
    int num_chromosomes = _genome.numChromosomes(0)+_genome.numChromosomes(1);
    color[] color_table = new color[num_chromosomes];      
    int counter = 0;                      
    for ( int  i = 0; i < _genome.numChromosomes(0); i++ ) {
       color_table[i] = color_map[i%color_map.length];
    }
    for ( int  i = 0; i < _genome.numChromosomes(1); i++ ) {
       color_table[_genome.genomeChromStartIndex(1)+i] = color_map[i%color_map.length];
    }
    
    // go through blocks and assign colors based on dest
    Block b;
    for ( int i = 0; i < num_chromosomes; i++ ) {
      for ( int j = 0; j < _genome.getChromosome(i).numSubunits(); j++ ) {
        b = _genome.getChromosome(i).getBlock(j);
        
        b.setColor( color_table[b.pair().refNumber()] );
        b.pair().setColor( color_table[b.pair().refNumber()] );
      }
    }
  } 
}
