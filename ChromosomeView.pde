float FRACTION_HEIGHT = 0.75;
float LETTER_LINE_H = 20.0;
int SCALE_OFFSET = 30;
int TRANSPOSON_COLOR = 0;
color _transposon_color = color( TRANSPOSON_COLOR );
float TRANSPOSON_LENGTH = 0.15;
float SELECTED_STROKE_WEIGHT = 2.0;
color SELECTED_COLOR = color( 0 );
int HISTO_BAR_WIDTH = 20;
int HISTO_BAR_OFFSET = 5;
int HISTO_BAR_HEIGHT = 4;

int _selected_block = 0;
boolean _block_selected;

float _max_block_pair_value;

boolean _block_rolled_over = false;
int _rolled_over_block = 0;

class ChromosomeView 
{
  Genome _genome;

  float _window_start, _window_end;
  int _width, _height;
  int _bar_width, _bar_height;
  float _scale;
  float _one_pixel_h;

  float _scale_unit;
  int _s_unit;

  int _text_box_x, _text_box_y, _text_box_w, _text_box_h;
  boolean _writing_text;
  String _input_text;

  ChromosomeView( Genome g, int w, int h ) {
    _genome = g;
    setDimensions( w, h );

    _window_start = 0.0;
    _window_end = 1.0;

    _block_selected = true;

    long min_len = _genome.getChromosome(_genome.genomeChromStartIndex(_selected_genome)).len();
    Float f_val;
    for ( int i = _genome.genomeChromStartIndex(_selected_genome)+1; i < _genome.genomeChromStartIndex(_selected_genome)+_genome.numChromosomes(_selected_genome); i++ ) 
    {
      f_val = min( min_len, _genome.getChromosome(i).len() );
      min_len = f_val.longValue();
    } 
    _scale_unit = min_len-(min_len%1.0e6);
    _s_unit = int( min_len / 1.0e6 );
    if ( _s_unit > 10 ) {
      _scale_unit = 10.0e6;
      _s_unit = 10;
    }
    else if ( _s_unit > 5 ) {
      _scale_unit = 5.0e6;
      _s_unit = 5;
    }
    else if ( _s_unit > 2 ) {
      _scale_unit = 2.0e6;
      _s_unit = 2;
    }
    setScale();

    // find the max variance of the blocks;
    _max_block_pair_value = Float.MIN_VALUE;
    for ( int i = 0; i < _genome.numChromosomes(0)+_genome.numChromosomes(1); i++ ) {
      for ( int j = 0; j < _genome.getChromosome(i).numSubunits(); j++ ) {       
        _max_block_pair_value = max( _max_block_pair_value, _genome.getChromosome(i).getBlock(j).pairValue() );
      }
    }

    _text_box_x = 0;
    _text_box_y = _bar_height + 20;
    _text_box_w = _bar_width;
    _text_box_h = 16;
    _writing_text = false;
    _input_text = "";
  }

  void setDimensions( int w, int h ) {
    _width = w;
    _height = h;

    _bar_width = w - HISTO_BAR_WIDTH - HISTO_BAR_OFFSET - SCALE_OFFSET;
    _bar_height = h;

    _one_pixel_h = 1.0/float(_bar_height);
  }

  void setScale() {
    _scale = _scale_unit/(float)_genome.getChromosome(_selected_chromosome).len();
  }

  void setWindowStartAndEnd( float start, float end ) {
    _window_start = max(min(start,1.0),0.0);
    _window_end = max(min(end,1.0),0.0);
  }

  void downKeyPressed() {
    if ( _block_selected ) {
      _selected_block = (_selected_block+1)%_genome.getChromosome(_selected_chromosome).numSubunits();
      _gene_selected = _gene_rolled_over = false;
      _selected_gene = _rolled_over_gene = 0;
    } 
  }

  void upKeyPressed() {
    if ( _block_selected ) {
      if ( _selected_block == 0 ) {
        _selected_block = _genome.getChromosome(_selected_chromosome).numSubunits() - 1;
      }
      else { 
        _selected_block = _selected_block - 1;
      }
      _gene_selected = _gene_rolled_over = false;
      _selected_gene = _rolled_over_gene = 0;
    } 
  }

  boolean overWindowShutterArea( int x, int y ) {
    return ( x <= SCALE_OFFSET );
  }

  boolean overStartShutter( float x, float y ) {
    float y_diff = abs((float)y - _window_start*_bar_height);

    if ( (x >= (SCALE_OFFSET-TRIANGLE_HEIGHT)) && (x <= SCALE_OFFSET) &&
      (y_diff <= HALF_TRIANGLE_HEIGHT) ) {
      return true; 
    }
    else {      
      return false;
    }
  }

  boolean overTextBox( float x, float y ) {
    if ( (x >= _text_box_x) && ( x <= _text_box_x+_text_box_w) &&
      (y >= _text_box_y) && ( y <= _text_box_y+_text_box_h) ) {
      _writing_text = true;
      _input_text = "";
      return true;
    }
    else {
      return false;
    }
  }

  void coordInput( char k ) {
    if ( _writing_text && (_input_text.length() < 10) ) {
      String s = Character.toString(k);
      _input_text = _input_text + s;
    }
  }

  void removeOneText() {
    if ( _writing_text && (_input_text.length() != 0 ) ) {
      _input_text = _input_text.substring( 0 , _input_text.length()-1 );
    }
  }

  void returnPressed() {
    _writing_text = false;

    // put commas in the input
    int val = parseInt( _input_text );
    if ( val != 0 ) {
      _input_text = nfc( val );
    }

    // find the closest block
    _selected_block = 0;
    long diff1, diff2;
    Float f_val;
    for ( int i = 1; i < _genome.getChromosome(_selected_chromosome).numSubunits(); i++ ) {
      if ( val <= _genome.getChromosome(_selected_chromosome).getBlock(i).stop() ) {
        diff1 = val - _genome.getChromosome(_selected_chromosome).getBlock(i-1).stop();
        f_val = abs(val - _genome.getChromosome(_selected_chromosome).getBlock(i).start());
        diff2 = f_val.longValue();

        if ( diff1 <= diff2 ) {
          _selected_block = i-1;
        }
        else {
          _selected_block = i;
        }
        return;
      }
    }
    _selected_block = _genome.getChromosome(_selected_chromosome).numSubunits() - 1;
  }

  boolean overEndShutter( float x, float y ) {
    float y_diff = abs((float)y - _window_end*_bar_height);

    if ( (x >= (SCALE_OFFSET-TRIANGLE_HEIGHT)) && (x <= SCALE_OFFSET) &&
      (y_diff <= HALF_TRIANGLE_HEIGHT) ) {
      return true; 
    }
    else {      
      return false;
    }
  }

  void updateShutter( float x, float y, boolean start_shutter, boolean end_shutter ) {
    float val = max(min((y / (float)_bar_height),1.0),0.0);   

    if ( start_shutter ) {
      _window_start = val;
    }
    else {
      _window_end = val;
    }
  }

  float windowStart() {
    return _window_start;
  }

  float windowEnd() {
    return _window_end;
  }

  boolean blockSelected( int x, int y ) {
    if ( (x >= SCALE_OFFSET) &&
      (x <= SCALE_OFFSET+_bar_width) )
    {
      // convert the y into a relative value
      float rel_y = (float)y/(float)_bar_height;

      // find out if this y value is within a region
      Block r;
      for ( int i = 0; i < _genome.getChromosome(_selected_chromosome).numSubunits(); i++ ) {
        r = _genome.getChromosome(_selected_chromosome).getBlock(i);
        if ( (rel_y >= r.relStart()) && (rel_y <= r.relStop()) ) {
          _selected_block = i;
          _block_selected = true;
          return true;
        }
      }
    }
    return false;
  }

  void setBlockSelected( boolean rs ) {
    _block_selected = rs;
  }

  boolean overBlock( int x, int y ) {
    if ( (x >= SCALE_OFFSET) &&
      (x <= SCALE_OFFSET+_bar_width) )
    {
      // convert the y into a relative value
      float rel_y = (float)y/(float)_bar_height;

      // find out if this y value is within a region
      Block r;
      for ( int i = 0; i < _genome.getChromosome(_selected_chromosome).numSubunits(); i++ ) {
        r = _genome.getChromosome(_selected_chromosome).getBlock(i);
        if ( (rel_y >= r.relStart()) && (rel_y <= r.relStop()) ) {
          _rolled_over_block = i;
          _block_rolled_over = true;
          return true;
        }
      }
    }
    _block_rolled_over = false;
    _rolled_over_block = 0;
    return false;
  }

  boolean isBlockSelected() {
    return _block_selected;
  }

  void render() {

    pushMatrix();
    translate(SCALE_OFFSET, 0 );

    // draw the window
    //noSmooth();
    noStroke(); 
    fill( 255 );    
    rectMode( CORNERS );
    boolean draw_all = false;
    boolean draw_partial = false;
    float s, e;
    if ( _window_start <= _window_end ){   
      s = min(max(_window_start, 0.0), 1.0);
      e = max(min(_window_end, 1.0),0.0);
      rect( 0, _bar_height*s, _bar_width, _bar_height*e );
    }
    else if ( (_window_start > _window_end) &&
      (((_window_start >= 1.0) && (_window_end >= 1.0)) ||
      ((_window_start <= 0.0) && (_window_end <= 0.0))) )
    {
      rect( 0, 0, _bar_width, _bar_height );
      draw_all = true;
    }
    else if ( (_window_start > _window_end) &&
      (_window_start >= 0.0) || (_window_end >=0.0) ||
      (_window_start <=1.0) || (_window_end <= 1.0 ) ) {
      e = max(_window_end, 0.0);
      rect( 0, 0, _bar_width, _bar_height*e );
      s = min(_window_start, 1.0);
      rect( 0, _bar_height*s, _bar_width, _bar_height );

      draw_partial = true;            
    }

    // draw the regions and blocks
    renderBlocks( draw_all, draw_partial );
    renderTracks( draw_all, draw_partial );

    // draw the tRNA
    //rendertRNA( draw_all, draw_partial );

    // draw the LTR
    //renderLTR( draw_all, draw_partial );

    // draw the line down the middle
    /*stroke( _background );
     strokeWeight( 1.0 );
     line( _bar_width/2.0, 0, _bar_width/2.0, _bar_height );*/

    // outline the block
    strokeWeight( 1 );
    stroke( 0 );
    noFill();
    rect( 0, 0, _bar_width, _bar_height );

    // draw the letters
    /*stroke( _chrom_color );
     line( _bar_width/2.0, -LETTER_LINE_H, _bar_width/2.0, 0 );
     float down = 5.0;
     smooth();
     fill( _chrom_color );
     textAlign( RIGHT, BOTTOM );
     text( "LTR", _bar_width/2.0 - down, -down );
     textAlign( LEFT, BOTTOM );
     text( "tRNA", _bar_width/2.0 + down, -down );*/

    // draw the scale bar
    //noSmooth();
    stroke( _chrom_color );
    strokeWeight( 2.0 );
    line( -SCALE_OFFSET, 0, -SCALE_OFFSET, _bar_height*_scale );
    line( -SCALE_OFFSET, 0, -SCALE_OFFSET*0.5, 0 );
    line( -SCALE_OFFSET, _bar_height*_scale, -SCALE_OFFSET*0.5, _bar_height*_scale );
    //smooth();
    fill( _chrom_color );
    textFont( sans_15 );
    textAlign( RIGHT, CENTER );
    text( _s_unit + "Mb", -SCALE_OFFSET-5.0, _bar_height*_scale );

    // if a region is selected, highlight it
    Block r;    
    if ( _block_rolled_over ) {
      r = _genome.getChromosome(_selected_chromosome).getBlock(_rolled_over_block);

      strokeWeight( SELECTED_STROKE_WEIGHT );
      stroke( SELECTED_COLOR );
      noFill();
      rectMode( CORNERS );
      rect( -1, _bar_height*r.relStart()-1, _bar_width+1, _bar_height*r.relStop()+1 );       
    }
    if ( _block_selected ) {
      r = _genome.getChromosome(_selected_chromosome).getBlock(_selected_block);

      strokeWeight( SELECTED_STROKE_WEIGHT );
      stroke( SELECTED_COLOR );
      noFill();
      rectMode( CORNERS );
      rect( -1, _bar_height*r.relStart()-1, _bar_width+1, _bar_height*r.relStop()+1 );    

      fill( SELECTED_COLOR );
      rect( -10, _bar_height*r.relStart()-1, -1, _bar_height*r.relStop()+1 );
    }

    renderNumber();
    renderSideHistogram();
    renderWindowSliders();

    renderInputBox();

    popMatrix();
  }

  void renderInputBox() {
    if ( _writing_text ) {
      fill( _button_color );
    }
    else {
      noFill();
    }
    stroke( 0 );
    rect( _text_box_x, _text_box_y, _text_box_w, _text_box_h );
    fill( 0 );
    textFont( sans_15 );
    textAlign( RIGHT, CENTER );
    text( "go to: ", _text_box_x, _text_box_y + _text_box_h/2 );

    textAlign( LEFT, CENTER );
    if ( _writing_text ) {
      fill( 0 );
    }
    else {
      fill( _chrom_color );
    }
    text( _input_text, _text_box_x+2, _text_box_y + _text_box_h/2 );
  }

  void renderWindowSliders() {
    // draw the triangle markers
    fill( color(0) );
    strokeWeight(1.0);
    stroke( _chrom_color );

    float y = _window_start*_bar_height;    
    triangle( 0, y, -TRIANGLE_HEIGHT, y-TRIANGLE_VERTEX_AXIS_OFFSET, -TRIANGLE_HEIGHT, y+TRIANGLE_VERTEX_AXIS_OFFSET );

    y = _window_end*_bar_height;    
    triangle( 0, y, -TRIANGLE_HEIGHT, y-TRIANGLE_VERTEX_AXIS_OFFSET, -TRIANGLE_HEIGHT, y+TRIANGLE_VERTEX_AXIS_OFFSET );
  }

  void renderSideHistogram() {
    //noSmooth();
    noStroke();
    fill( 255 );
    rectMode( CORNER );
    rect( _bar_width+HISTO_BAR_OFFSET, 0, HISTO_BAR_WIDTH, _bar_height );

    Block b;
    int half_histo_height = HISTO_BAR_HEIGHT/2;
    float half_way;
    fill( 0 );
    noStroke();
    for ( int i = 0; i < _genome.getChromosome(_selected_chromosome).numSubunits(); i++ ) {  
      b = _genome.getChromosome(_selected_chromosome).getBlock(i);     
      half_way = b.relStart() + (b.relStop() - b.relStart())/2.0;

      rect( _bar_width+HISTO_BAR_OFFSET, half_way*_bar_height-half_histo_height, b.pairValue()/_max_block_pair_value*HISTO_BAR_WIDTH, HISTO_BAR_HEIGHT );
    }

    /*float x = _bar_width+HISTO_BAR_OFFSET;
     float y = _bar_height;
     fill( _chrom_color );
     textAlign( CENTER, TOP );
     text( "similarity", (x+HISTO_BAR_WIDTH/2), (y+5) ); 
     text( "variance", (x+HISTO_BAR_WIDTH/2), (y+18) ); */
  }

  void renderNumber() {
    //smooth();
    noStroke();
    textFont( sans_12 );
    textAlign( CENTER, CENTER );

    float s = 0.5*_bar_width;
    float h = -40;

    fill( 255 );
    ellipse( s, h, TEXT_DIAMETER, TEXT_DIAMETER );
    fill( 0 );   
    text( _genome.getChromosome(_selected_chromosome).tag(), s, h );
  }

  void renderBlocks( boolean draw_all, boolean draw_partial ) {
    Block r;
    boolean use_region_color;
    //rectMode( CORNER );
    float s;
    for ( int i = 0; i < _genome.getChromosome(_selected_chromosome).numSubunits(); i++ ) {  
      r = _genome.getChromosome(_selected_chromosome).getBlock(i); 

      if ( !draw_all && !draw_partial ) {
        if ( (r.relStart() < _window_start) || (r.relStop() > _window_end) ) {
          use_region_color = false;          
        }
        else {
          use_region_color = true;
        }
      }
      else if ( draw_partial ) {
        if ( (r.relStop() > _window_end) && (r.relStart() < _window_start) ) {
          use_region_color = false;
        }
        else {
          use_region_color = true;
        }
      }  
      else {
        use_region_color = true;
      } 

      if ( use_region_color ) {
        fill( r.getColor() );
      }
      else {
        fill( _chrom_color );
      }

      rectMode( CORNER );
      rect( 0, r.relStart()*_bar_height, _bar_width, max(_bar_height*(r.relStop()-r.relStart()),1.0) ); 
    }  
  }

  void renderTracks( boolean draw_all, boolean draw_partial ) {
    int interval = _bar_width / (_genome.numTracks()+1);
    int x = interval;
    Annotation annot;
    float start;
    strokeWeight( 4 );
    strokeCap(SQUARE);
    fill( 0 );
    textFont( sans_12 );
    textAlign( CENTER, BOTTOM );
    for ( int t = 0; t < _genome.numTracks(); t++ ) {
      for ( int a = 0; a < _genome.getChromosome(_selected_chromosome).numAnnotations(t); a++ ) {
        annot = _genome.getChromosome(_selected_chromosome).getAnnotation(t,a);

        if ( !draw_all && !draw_partial ) {
          if ( (annot.relStart() < _window_start) || (annot.relStop() > _window_end) ) {
            stroke( _chrom_color );          
          }
          else {
            stroke( 0 );
          }
        }
        else if ( draw_partial ) {
          if ( (annot.relStop() > _window_end) && (annot.relStart() < _window_start) ) {
            stroke( _chrom_color ); 
          }
          else {
            stroke( 0 );
          }
        }  
        else {
          stroke( 0 );
        } 
        
        start = annot.relStart()*_bar_height;
        line( x, start, x, start+max(_bar_height*(annot.relStop()-annot.relStart()),1.0) );
      }      
      text( _genome.trackTag(t), x, -5 );
      
      x += interval;
    }
    
  }

  /*  void rendertRNA( boolean draw_all, boolean draw_partial ) {
   Transposon t;
   boolean use_transposon_color;
   strokeWeight( 1.0 );
   int tw = int(_bar_width*TRANSPOSON_LENGTH);
   for ( int i = 0; i < _genome.getChromosome(_selected_chromosome).numtRNA(); i++ ) {
   t = _genome.getChromosome(_selected_chromosome).gettRNA(i); 
   
   if ( !draw_all && !draw_partial ) {
   if ( (t.relStart() < _window_start) || (t.relStop() > _window_end) ) {
   use_transposon_color = false;          
   }
   else {
   use_transposon_color = true;
   }
   }
   else if ( draw_partial ) {
   if ( (t.relStop() > _window_end) && (t.relStart() < _window_start) ) {
   use_transposon_color = false;
   }
   else {
   use_transposon_color = true;
   }
   }  
   else {
   use_transposon_color = true;
   }
   
   if ( use_transposon_color ) {
   stroke( _transposon_color );
   }
   else {
   stroke( _chrom_color );
   }
   line( _bar_width/2.0, t.relStart()*_bar_height, (_bar_width/2.0 + tw), t.relStart()*_bar_height );       
   }
   }
   
   void renderLTR( boolean draw_all, boolean draw_partial ) {
   Transposon t;
   boolean use_transposon_color;
   strokeWeight( 1.0 );
   int tw = int(_bar_width*TRANSPOSON_LENGTH);
   for ( int i = 0; i < _genome.getChromosome(_selected_chromosome).numLTR(); i++ ) {
   t = _genome.getChromosome(_selected_chromosome).getLTR(i); 
   
   if ( !draw_all && !draw_partial ) {
   if ( (t.relStart() < _window_start) || (t.relStop() > _window_end) ) {
   use_transposon_color = false;          
   }
   else {
   use_transposon_color = true;
   }
   }
   else if ( draw_partial ) {
   if ( (t.relStop() > _window_end) && (t.relStart() < _window_start) ) {
   use_transposon_color = false;
   }
   else {
   use_transposon_color = true;
   }
   }  
   else {
   use_transposon_color = true;
   }
   
   if ( use_transposon_color ) {
   stroke( _transposon_color );
   }
   else {
   stroke( _chrom_color );
   }
   line( _bar_width/2.0, t.relStart()*_bar_height, (_bar_width/2.0 - tw), t.relStart()*_bar_height );       
   }
   }*/
}








