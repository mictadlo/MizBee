
int ARROW_HEAD_LENGTH = 5;
int SPACE_BETWEEN = 200;
color _positive_orientation_color = color(220, 180, 180 );
color _negative_orientation_color = color( 200, 200, 220 );

color _button_color = color( 137, 205, 255 );

boolean _gene_selected;
int _selected_gene;

boolean _gene_rolled_over;
int _rolled_over_gene;

float _max_pair_value;

class BlockView 
{
  int _width, _height, _bar_width, _bar_height, _space_inbetween;
  Genome _genome;
  int m_x, m_y;
  boolean _over_selected_block;
  int _histo_l, _histo_offset;
  float _max_pair_value;
  
  boolean _invert;
  int _invert_x, _invert_y, _invert_w;
  
  int _num_zoom_levels, _selected_zoom_level;
  float[] _zoom_sizes;
  boolean _have_zoom;
  float _zoom_center, _zoom_control_step;
  int _zoom_bar_y, _zoom_bar_r, _zoom_bar_w;
  
  int _scroll_bar_x, _scroll_bar_w;
  float _scroll_bar_half_h;

  BlockView( Genome g, int w, int h ) {
    _genome = g;
    setDimensions( w, h );

    _gene_rolled_over = false;
    _rolled_over_gene = 0;
    _gene_selected = false;
    _selected_gene = 0;
    m_x = m_y = 0;
    
    _invert = false;
    
    

    // find the max variance of the blocks;
    _max_pair_value = Float.MIN_VALUE;
    for ( int i = 0; i < _genome.numChromosomes(0)+_genome.numChromosomes(1); i++ ) {
      for ( int k = 0; k < _genome.getChromosome(i).numSubunits(); k++ ) {
        for ( int l = 0; l < _genome.getChromosome(i).getBlock(k).numSubunits(); l++ ) {
          _max_pair_value = max( _max_pair_value, _genome.getChromosome(i).getBlock(k).getGene(l).pairValue() );
        }
      }
    }
  }

  void setDimensions( int w, int h ) {
    _width = w;
    _height = h;
  
    _space_inbetween = min( SPACE_BETWEEN, int( (float)((w-2*(HISTO_BAR_WIDTH+HISTO_BAR_OFFSET)))*0.65 ) );

    _scroll_bar_w = HISTO_BAR_WIDTH/2;
    int sb_offset = HISTO_BAR_OFFSET/4;
    _scroll_bar_x = w - _scroll_bar_w - HISTO_BAR_WIDTH+HISTO_BAR_OFFSET;
    
    _histo_l = 3*_width/4;    
    _bar_width = (_width - _space_inbetween - 2*(HISTO_BAR_WIDTH+HISTO_BAR_OFFSET) - _scroll_bar_w - sb_offset)/2;   
    _invert_w = _bar_width;
    _zoom_bar_r = 10;
    
    _zoom_bar_w = 2*_bar_width + _space_inbetween;
    _invert_x = _bar_width+_space_inbetween;
    
    int inbetween = 20;
       
    _bar_height = h-_histo_l-_invert_w - 4*inbetween;    
    _zoom_bar_y = _bar_height + 2*inbetween;
    _invert_y = _zoom_bar_y + inbetween;
    _histo_offset = _invert_w + 4*inbetween;

    //setUpZoom( _genome.minGeneBlockRatio() );  
    setUpZoom( _genome.getChromosome(_selected_chromosome).getBlock(_selected_block).minGeneBlockRatio() );
  }
  
  void setUpZoom( float min_val ) {
    // create the zoom levels
    float rel_length_per_pixel = 1.0/(float)_bar_height;
    float min_num_pixels = 5.0;
    
    // find out how many pixels we have for the smallest gene
    float pix_per_min_gene = min_val / rel_length_per_pixel;

    // if the smallest gene is big enough, don't enable zooming
    if ( pix_per_min_gene >= min_num_pixels ) {
      _num_zoom_levels = 1;
      _have_zoom = false;      
    }
    else {
      _num_zoom_levels = 4;
      _have_zoom = true;
    }
    
    // for now, do a constant amount of zoom wrt the total block size
    float zoom_step_size;
    if ( _num_zoom_levels > 1 ) {
      zoom_step_size = (1.0 - (pix_per_min_gene/min_num_pixels)) / float(_num_zoom_levels-1);
    }
    else {
      zoom_step_size = 0.0;
    }
    _zoom_sizes = new float[_num_zoom_levels];
    for ( int i = 0; i < _num_zoom_levels; i++ ) {
      _zoom_sizes[i] = 1.0 - (float)i*zoom_step_size;
    } 
    
    _selected_zoom_level = 0;  
    _zoom_center = 0.5;  
    
    _zoom_control_step = (float)_zoom_bar_w / float(_num_zoom_levels-1);
    
    setScrollBarHalfH();
  }
  
  void setScrollBarHalfH() {
    _scroll_bar_half_h = _zoom_sizes[_selected_zoom_level]*_bar_height/2;
  }

  boolean overGene( int x, int y ) {
    // convert the y into a relative value
    float rel_y = (float)y/(float)_bar_height;
    
    if ( (rel_y < 0.0) || (rel_y > 1.0) ) {
      _gene_selected = false;
      _gene_rolled_over = false;
      return false;
    }
    
    // convert to the relative global space
    rel_y *= _zoom_sizes[_selected_zoom_level];
    
    // account for center
    rel_y += (_zoom_center - _zoom_sizes[_selected_zoom_level]/2.0);
    
    
    //float rel_y = (float)y/(float)_bar_height;
    

    // determine which block we are over
    if ( (x >= HISTO_BAR_OFFSET+HISTO_BAR_WIDTH) && (x <= HISTO_BAR_OFFSET+HISTO_BAR_WIDTH+_bar_width) ) {
      _over_selected_block = true;
    }
    else if ( (x >= (_bar_width + _space_inbetween +HISTO_BAR_OFFSET+HISTO_BAR_WIDTH)) && 
              (x <= 2.0*_bar_width + _space_inbetween + HISTO_BAR_OFFSET+HISTO_BAR_WIDTH) ) {
      _over_selected_block = false;
      if ( _invert ) {
        rel_y = 1.0 - rel_y;
      }
    }
    else {
      _gene_rolled_over = false;
      _gene_selected = false;
      return false;
    }

    // find out if this y value is within a gene
    Gene g;
    for ( int i = 0; i < _genome.getChromosome(_selected_chromosome).getBlock(_selected_block).numSubunits(); i++ ) {
      if ( _over_selected_block ) {
        g = _genome.getChromosome(_selected_chromosome).getBlock(_selected_block).getGene(i);
      }
      else {
        g = _genome.getChromosome(_selected_chromosome).getBlock(_selected_block).getGene(i).pair();
        
      }
      
      if ( (rel_y >= g.containerRelStart()) && (rel_y <= g.containerRelStop()) ) {
        _rolled_over_gene = i;
        _gene_rolled_over = true;
        m_x = x - HISTO_BAR_OFFSET-HISTO_BAR_WIDTH; 
        m_y = y;
        _gene_selected = true;
        _selected_gene = i;
        return true;
      }
    }
    _gene_rolled_over = false;
    _gene_selected = false;
    return false;
  }
  
  boolean invertSelected( int x, int y ) {
    x -= HISTO_BAR_OFFSET+HISTO_BAR_WIDTH;
    if ( (x >= _invert_x) && (x <= _invert_x+_invert_w) && (y >= _invert_y) && (y <= _invert_y+_invert_w) ) {
      _invert = !_invert;
      return true;
    }
    return false;
  }
  
  boolean overZoomControl( int x, int y ) {
    if ( !_have_zoom ) {
      return false;
    }
    
    x -= HISTO_BAR_OFFSET+HISTO_BAR_WIDTH;
    if ( (x >= 0) && (x <= _zoom_bar_w) && (y >= _zoom_bar_y-_zoom_bar_r) && ( y <= _zoom_bar_y+ _zoom_bar_r) ) {
      _selected_zoom_level = round( (float)x / _zoom_control_step );  
      setScrollBarHalfH();
      scrollTheBar( int(_zoom_center*_bar_height) );
      return true;
    }    
    return false;
  }
  
  void zoom() {
    _selected_zoom_level = min( ++_selected_zoom_level, _num_zoom_levels-1 );
    setScrollBarHalfH();
  }
    
  boolean overScrollControl( int x, int y ) { 
    x -= HISTO_BAR_WIDTH+HISTO_BAR_OFFSET;
    if ( (x >= _scroll_bar_x) && (x <= _scroll_bar_x+_scroll_bar_w) ) {
      float center = _zoom_center*_bar_height;
      if ( (y >= center-_scroll_bar_half_h) && (y <= center+_scroll_bar_half_h) ) {
        return true;
      }
    }
    return false;
  }
  
  void scroll( int y ) {
    //y = 0;
    // but y in space relative to the size the frame [0,1]
    float y_val = (float)y/(float)_bar_height;
    
    // convert to the relative global space
    y_val *= _zoom_sizes[_selected_zoom_level];
    
    // account for center
    y_val += (_zoom_center - _zoom_sizes[_selected_zoom_level]/2.0);
    
    // now, back to pixel space
    y_val *= (float)_bar_height;
    
    _selected_zoom_level = min( ++_selected_zoom_level, _num_zoom_levels-1 );
    setScrollBarHalfH();

    float center = min( max(y_val, _scroll_bar_half_h), _bar_height - _scroll_bar_half_h );
    _zoom_center = center / (float)_bar_height;
  }
  
  void scrollTheBar( int y ) {
    float center = min( max((float)y, _scroll_bar_half_h), _bar_height - _scroll_bar_half_h );
    _zoom_center = center / (float)_bar_height;
  }
    

  boolean geneSelected( int x, int y ) {
    // convert the y into a relative value
    float rel_y = (float)y/(float)_bar_height;
    
    if ( (rel_y < 0.0) || (rel_y > 1.0) ) {
      _gene_selected = false;
      return false;
    }
    
    // convert to the relative global space
    rel_y *= _zoom_sizes[_selected_zoom_level];
    
    // account for center
    rel_y += (_zoom_center - _zoom_sizes[_selected_zoom_level]/2.0);
    
    
    //float rel_y = (float)y/(float)_bar_height;
    

    // determine which region we are over
    boolean over_selected_block;
    if ( (x >= HISTO_BAR_OFFSET+HISTO_BAR_WIDTH) && (x <= _bar_width + HISTO_BAR_OFFSET+HISTO_BAR_WIDTH) ) {
      over_selected_block = true;
    }
    else if ( (x >= (_bar_width + _space_inbetween + HISTO_BAR_OFFSET+HISTO_BAR_WIDTH)) && 
              (x <= 2.0*_bar_width + _space_inbetween + HISTO_BAR_OFFSET+HISTO_BAR_WIDTH) ) {
      over_selected_block = false;
      if ( _invert ) {
        rel_y = 1.0 - rel_y;
      }
    }
    else {
      _gene_selected = false;
      return false;
    }

    // find out if this y value is within a region
    Gene g;
    for ( int i = 0; i < _genome.getChromosome(_selected_chromosome).getBlock(_selected_block).numSubunits(); i++ ) {
      if ( over_selected_block ) {
        g = _genome.getChromosome(_selected_chromosome).getBlock(_selected_block).getGene(i);
      }
      else {
        g = _genome.getChromosome(_selected_chromosome).getBlock(_selected_block).getGene(i).pair();
        
      }

      if ( (rel_y >= g.containerRelStart()) && (rel_y <= g.containerRelStop()) ) {
        _selected_gene = i;
        _gene_selected = true;
        return true;
      }
    }
    _gene_selected = false;
    return false;
  }

  boolean isGeneSelected() {
    return _gene_selected;
  }

  void setGeneSelected( boolean s ) {
    _gene_selected = s;
  }

  void render() {
    pushMatrix();
    translate( HISTO_BAR_OFFSET+HISTO_BAR_WIDTH, 0 );
    
    if ( _block_selected ) {
      //noSmooth();
      rectMode( CORNER );

      Block b, p;
      b = _genome.getChromosome(_selected_chromosome).getBlock(_selected_block);
      p = _genome.getChromosome(_selected_chromosome).getBlock(_selected_block).pair();

      // render the underlying region
      /*noStroke();
       fill( b.getColor() );
       rect( 0, 0, _bar_width, _height );
       fill( p.getColor() );
       rect( _bar_width+_space_inbetween, 0, _bar_width, _height );*/

      // render the outline           
      //noFill();     
      fill( lerpColor( b.getColor(), #FFFFFF, 0.8 ) );
      stroke( _chrom_color );
      strokeWeight( 1 );
      rect( 0, 0, _bar_width, _bar_height );
      rect( _bar_width+_space_inbetween, 0, _bar_width, _bar_height );

      

      // render the blocks in this region and its pair
      float s, h, ps, ph;
      float x, y, h_r, h_g;
      Long l_val;
      for ( int i = 0; i < b.numSubunits(); i++ ) {
        l_val = b.getGene(i).start()-b.start();
        s = l_val.floatValue() / (float)b.len();
        h = (float)b.getGene(i).len() / (float)b.len();

        l_val = p.getGene(i).start()-p.start();
        ps = l_val.floatValue() / (float)p.len();
        ph = (float)p.getGene(i).len() / (float)p.len();
        
        if ( _invert ) {
          ps = 1.0 - ps;
        }
        
        s -= (_zoom_center-_zoom_sizes[_selected_zoom_level]/2.0);
        s /= _zoom_sizes[_selected_zoom_level];        
        h /= _zoom_sizes[_selected_zoom_level];
        ps -= (_zoom_center-_zoom_sizes[_selected_zoom_level]/2.0);
        ps /= _zoom_sizes[_selected_zoom_level];        
        ph /= _zoom_sizes[_selected_zoom_level];

        // draw the connecting lines 
        //smooth();    
        if ( _gene_selected && (_selected_gene == i) ) {
          strokeWeight( 1 );
          stroke( SELECTED_COLOR );
        }
        else {  
          strokeWeight( 1 );
          stroke( color(_background) );
        }
        if ( b.getGene(i).pairOrientation() ) {
          fill( _positive_orientation_color );         
        }
        else {
          fill( _negative_orientation_color );
        }
        if ( !_invert ) {
          quad( _bar_width-1, s*_bar_height-1, _bar_width+_space_inbetween, ps*_bar_height-1, 
          _bar_width+_space_inbetween, (ps+ph)*_bar_height+1, _bar_width-1, (s+h)*_bar_height+1 );
        }
        else {
          quad( _bar_width-1, s*_bar_height-1, _bar_width+_space_inbetween, (ps-ph)*_bar_height+1, _bar_width+_space_inbetween, ps*_bar_height+1, 
          _bar_width-1, (s+h)*_bar_height+1 );
        }
        
        strokeWeight( 1.0 );
        if ( b.getGene(i).pairOrientation() ) {
          stroke( _positive_orientation_color );         
        }
        else {
          stroke( _negative_orientation_color );
        }
        if ( !_invert ) {
          line( _bar_width, (s+h)*_bar_height, _bar_width+_space_inbetween, (ps+ph)*_bar_height );
        }
        else {
          line( _bar_width, (s+h)*_bar_height, _bar_width+_space_inbetween, (ps-ph)*_bar_height );
        }

        // draw the gene
        //noSmooth();
        //smooth();

        if ( _gene_selected && (_selected_gene == i) ) {
          strokeWeight( SELECTED_STROKE_WEIGHT );
          stroke( SELECTED_COLOR );
        }
        else {  
          strokeWeight( 1 );
          stroke( 0 );
        }         

        fill( b.getColor() );

        if ( b.getGene(i).orientation() ) {
          x = 0;
          y = s*_bar_height;
          h_g = h*_bar_height;
          h_r = h_g - ARROW_HEAD_LENGTH;
          if ( h_r > 0 ) {
            beginShape();
            vertex( x, y );
            vertex( x+_bar_width, y );
            vertex( x+_bar_width, y+h_r );
            vertex( x+_bar_width/2, y+h_g );
            vertex( x, y+h_r );
            endShape(CLOSE);
          }
          else {
            beginShape();
            vertex( x, y );
            vertex( x+_bar_width, y );
            vertex( x+_bar_width/2, y+h_g );
            endShape(CLOSE);
          }
        }
        else {
          x = 0;
          y = s*_bar_height;
          h_g = h*_bar_height;
          h_r = ARROW_HEAD_LENGTH;
          if ( h_r < h_g ) {
            beginShape();
            vertex( x+_bar_width/2, y );
            vertex( x+_bar_width, y+h_r );
            vertex( x+_bar_width, y+h_g );
            vertex( x, y+h_g );
            vertex( x, y+h_r );
            endShape(CLOSE);
          }
          else {
            beginShape();
            vertex( x+_bar_width/2, y );
            vertex( x+_bar_width, y+h_g );
            vertex( x, y+h_g );
            endShape(CLOSE);
          }
        }

        fill( b.getColor() );

        if ( b.getGene(i).pair().orientation() ) {
          x = _bar_width+_space_inbetween;
          y = ps*_bar_height;
          h_g = ph*_bar_height;
          h_r = h_g - ARROW_HEAD_LENGTH;
          if ( !_invert ) {
            if ( h_r > 0 ) {
              beginShape();
              vertex( x, y );
              vertex( x+_bar_width, y );
              vertex( x+_bar_width, y+h_r );
              vertex( x+_bar_width/2, y+h_g );
              vertex( x, y+h_r );
              endShape(CLOSE);
            }
            else {
              beginShape();
              vertex( x, y );
              vertex( x+_bar_width, y );
              vertex( x+_bar_width/2, y+h_g );
              endShape(CLOSE);
            }
          }
          else {
            if ( h_r > 0 ) {
              beginShape();
              vertex( x, y );
              vertex( x+_bar_width, y );
              vertex( x+_bar_width, y-h_r );
              vertex( x+_bar_width/2, y-h_g );
              vertex( x, y-h_r );
              endShape(CLOSE);
            }
            else {
              beginShape();
              vertex( x, y );
              vertex( x+_bar_width, y );
              vertex( x+_bar_width/2, y-h_g );
              endShape(CLOSE);
            }
          }
        }
        else {
          x = _bar_width+_space_inbetween;
          y = ps*_bar_height;
          h_g = ph*_bar_height;
          h_r = ARROW_HEAD_LENGTH;
          if ( !_invert ) {
            if ( h_r < h_g ) {
              beginShape();
              vertex( x+_bar_width/2, y );
              vertex( x+_bar_width, y+h_r );
              vertex( x+_bar_width, y+h_g );
              vertex( x, y+h_g );
              vertex( x, y+h_r );
              endShape(CLOSE);
            }
            else {
              beginShape();
              vertex( x+_bar_width/2, y );
              vertex( x+_bar_width, y+h_g );
              vertex( x, y+h_g );
              endShape(CLOSE);
            }
          }
          else {
            if ( h_r < h_g ) {
              beginShape();
              vertex( x+_bar_width/2, y );
              vertex( x+_bar_width, y-h_r );
              vertex( x+_bar_width, y-h_g );
              vertex( x, y-h_g );
              vertex( x, y-h_r );
              endShape(CLOSE);
            }
            else {
              beginShape();
              vertex( x+_bar_width/2, y );
              vertex( x+_bar_width, y-h_g );
              vertex( x, y-h_g );
              endShape(CLOSE);
            }
          }
        }
      }
      
      renderSideHistograms();
      
      // clean up from the zoom
      fill( _background );
      noStroke();
      rect( -2-HISTO_BAR_OFFSET-HISTO_BAR_WIDTH, 0, _w+2, -1000 );
      rect( -2-HISTO_BAR_OFFSET-HISTO_BAR_WIDTH, _bar_height, _w+2, 1000 );
      
      renderHistogram();
      renderNumbers( _genome.getChromosome(b.refNumber()).tag(), _genome.getChromosome(p.refNumber()).tag() );

      // render the label
      /*smooth();
      noStroke();
      textAlign( CENTER, BOTTOM );
      fill( _chrom_color );
      text( "block", _width/2-HISTO_BAR_OFFSET-HISTO_BAR_WIDTH, -70 );*/

      if ( _gene_rolled_over ) {
        renderGeneID();
      }
      
      renderOrientationLegend();
      renderInvert();
      
      if ( _have_zoom ) {
        renderScrollBar();
        renderZoomBar();
      }
    }
    
    popMatrix();
  }
  
  void renderScrollBar() {
    // render the bar
    noFill();
    stroke( _chrom_color );
    rect( _scroll_bar_x, 0, _scroll_bar_w, _bar_height );
    
    fill( _button_color );
    stroke( 0 );
    rectMode( CENTER );
    rect( _scroll_bar_x+_scroll_bar_w/2, _bar_height*_zoom_center, _scroll_bar_w, _bar_height*_zoom_sizes[_selected_zoom_level] );
    rectMode( CORNER );
  }
  
  void renderZoomBar() {
    stroke( 0 );
    line( 0, _zoom_bar_y, _zoom_bar_w, _zoom_bar_y );
    line( 0, _zoom_bar_y - _zoom_bar_r, 0, _zoom_bar_y + _zoom_bar_r );
    line( _zoom_bar_w, _zoom_bar_y - _zoom_bar_r, _zoom_bar_w, _zoom_bar_y + _zoom_bar_r );
    
    
    int half_r = _zoom_bar_r/2;
    for ( int i = 1; i < _num_zoom_levels; i++ ) {
      line( _zoom_control_step*i, _zoom_bar_y-half_r, _zoom_control_step*i, _zoom_bar_y+half_r );
    }
    
    fill( _button_color );
    stroke( 0 );
    ellipse( _zoom_control_step*_selected_zoom_level, _zoom_bar_y, _zoom_bar_r, _zoom_bar_r );
    
    textFont( sans_15 );
    textAlign( RIGHT, CENTER );
    fill( 0 );
    text( "out", -10, _zoom_bar_y );
    textAlign( LEFT, CENTER );
    text( "in", _zoom_bar_w+10, _zoom_bar_y );
  }
  
  void renderInvert() {
    
    if ( _invert ) {
      fill( _button_color );
    }
    else {
      noFill();
    }
    strokeWeight( 1 );
    //stroke( _chrom_color );
    stroke( 0 );
    rect( _invert_x, _invert_y, _invert_w, _invert_w );
    
    textFont( sans_15 );
    textAlign( LEFT, CENTER );
    //fill( _chrom_color );
    fill( 0 );
    text( " invert", _invert_x+_invert_w, _invert_y+_invert_w/2 );
  }

  void renderSideHistograms() { 
    float x = 2.0*_bar_width+_space_inbetween;

    //noSmooth();
    noStroke();
    fill( 255 );
    rectMode( CORNER );
    rect( -HISTO_BAR_OFFSET-HISTO_BAR_WIDTH, 0, HISTO_BAR_WIDTH, _bar_height );
    rect( x+HISTO_BAR_OFFSET, 0, HISTO_BAR_WIDTH, _bar_height );

    Block b, p;
    Gene g;
    int half_histo_height = HISTO_BAR_HEIGHT/2;
    float half_way;
    fill( 0 );
    noStroke();

    b = _genome.getChromosome(_selected_chromosome).getBlock(_selected_block);     
    p = b.pair();
    float s, h;
    Long l_val;
    for ( int j = 0; j < b.numSubunits(); j++ ) {
      g = b.getGene(j);

      l_val = g.start()-b.start();
      s = l_val.floatValue() / (float)b.len();
      h = (float)g.len() / (float)b.len();     
      half_way = s + h/2.0;    
    
      s -= (_zoom_center-_zoom_sizes[_selected_zoom_level]/2.0);
      s /= _zoom_sizes[_selected_zoom_level];        
      h /= _zoom_sizes[_selected_zoom_level]; 
      half_way -= (_zoom_center-_zoom_sizes[_selected_zoom_level]/2.0);
      half_way /= _zoom_sizes[_selected_zoom_level]; 
      
      rect( -HISTO_BAR_OFFSET, half_way*_bar_height-half_histo_height, -(g.pairValue()/_max_pair_value*HISTO_BAR_WIDTH), HISTO_BAR_HEIGHT );

      l_val = g.pair().start()-p.start();
      s = l_val.floatValue() / (float)p.len();
      h = (float)g.pair().len() / (float)p.len();
      if ( _invert ) {
        s = 1.0 - s;
        half_way = s - h/2.0;
      }
      else {
        half_way = s + h/2.0;
      }      
      
      s -= (_zoom_center-_zoom_sizes[_selected_zoom_level]/2.0);
      s /= _zoom_sizes[_selected_zoom_level];        
      h /= _zoom_sizes[_selected_zoom_level]; 
      half_way -= (_zoom_center-_zoom_sizes[_selected_zoom_level]/2.0);
      half_way /= _zoom_sizes[_selected_zoom_level]; 
      
      rect( x+HISTO_BAR_OFFSET, half_way*_bar_height-half_histo_height, g.pairValue()/_max_pair_value*HISTO_BAR_WIDTH, HISTO_BAR_HEIGHT );
    }
  }

  void renderSideHistogram() { 
    float x = 2.0*_bar_width+_space_inbetween;

    //noSmooth();
    noStroke();
    fill( 255 );
    rectMode( CORNER );
    rect( x+HISTO_BAR_OFFSET, 0, HISTO_BAR_WIDTH, _bar_height );

    Block b, p;
    Gene g;
    int half_histo_height = HISTO_BAR_HEIGHT/2;
    float half_way;
    fill( 0 );
    noStroke();

    b = _genome.getChromosome(_selected_chromosome).getBlock(_selected_block);     
    p = b.pair();
    float s, h;
    Long l_val;
    for ( int j = 0; j < b.numSubunits(); j++ ) {
      g = b.getGene(j);
      l_val = g.pair().start()-p.start();
      s = l_val.floatValue() / (float)p.len();
      h = (float)g.pair().len() / (float)p.len();
      half_way = s + h/2.0;

      rect( x+HISTO_BAR_OFFSET, half_way*_bar_height-half_histo_height, g.pairValue()/_max_pair_value*HISTO_BAR_WIDTH, HISTO_BAR_HEIGHT );
    }
  }

  void renderHistogram() {   
    float x = float( (_width-_histo_l)/2 )-HISTO_BAR_OFFSET-HISTO_BAR_WIDTH;
    float y = float(_bar_height+_histo_offset);

    fill(255);
    noStroke();
    rectMode( CORNER );
    rect( x, y, _histo_l, _histo_l ); 

    Block b = _genome.getChromosome(_selected_chromosome).getBlock(_selected_block);
    int num_bars = b.numSubunits();
    float histo_bar_w = (float)_histo_l / (float)num_bars;
    for ( int i = 0; i < num_bars; i++ ) {
      if ( _gene_rolled_over && (_rolled_over_gene == i) ) {
        textFont( sans_12 );
        textAlign( RIGHT, CENTER );
        fill( 0 );
        String s = nf( b.getGene(i).pairValue(), 1, 2 );
        text( s, (x-5), (y+0.5*histo_bar_w) );
        noStroke();
        //fill( 0 );
      }
      else {
        //stroke( _background );
        fill( _chrom_color );
      }
      rect( x, y, (b.getGene(i).pairValue()/_max_pair_value)*_histo_l, histo_bar_w );
      y += histo_bar_w;
    }

    /*fill( _chrom_color );
    textAlign( CENTER, TOP );
    text( "similarity", (x+_histo_l/2), (y+5) ); 
    text( "score", (x+_histo_l/2), (y+18) ); */
  }

  void renderGeneID() {
    String id1;
    if ( _over_selected_block ) {
      id1 = _genome.getChromosome(_selected_chromosome).getBlock(_selected_block).getGene(_rolled_over_gene).tag();
    }
    else {
      id1 = _genome.getChromosome(_selected_chromosome).getBlock(_selected_block).getGene(_rolled_over_gene).pair().tag();
    }
    if ( id1.length() == 0 ) {
      return;
    }
      
    //smooth();
    noStroke();
    textFont( sans_12 );
    textAlign( CENTER, CENTER );

    int w = id1.length()*10;
    int h = 20;
    int x = m_x /*+ w/2*/;
    int y = m_y - h/2 - 5;

    fill( 255 );
    rectMode( CENTER );
    rect( x, y, w, h );
    fill( 0 );
    if ( _over_selected_block ) {
      text( id1, x, y );
    }
    else {
      text( id1, x, y );
    }    
  }

  void renderNumbers( String r_num, String p_num ) {
    //smooth();
    noStroke();
    textFont( sans_12);
    textAlign( CENTER, CENTER );

    float s = 0.5*_bar_width;
    float ps = 1.5*_bar_width + _space_inbetween;
    float h = -40;

    fill( 255 );
    ellipse( s, h, TEXT_DIAMETER, TEXT_DIAMETER );
    //ellipse( ps, h, TEXT_DIAMETER, TEXT_DIAMETER );

    fill( 0 );   
    text( r_num, s, h );
    text( p_num, ps, h );
    
    float sp_start = _zoom_center - _zoom_sizes[_selected_zoom_level]/2;
    float sp_stop = _zoom_center + _zoom_sizes[_selected_zoom_level]/2;
    Long coord;
    
    
    h = -5;
    int text_offset = _space_inbetween/2;
    fill( 0 );
    Long l_val;
    Float f_val;
    l_val = _genome.getChromosome(_selected_chromosome).getBlock(_selected_block).len();
    f_val = l_val.floatValue()*sp_start;
    coord = _genome.getChromosome(_selected_chromosome).getBlock(_selected_block).start()+f_val.longValue();
    //coord -= int(_zoom_center-_zoom_sizes[_selected_zoom_level]/2.0);
    //coord /= _zoom_sizes[_selected_zoom_level];
    textAlign( RIGHT, BOTTOM );
    text( coord.toString(), s+text_offset, h );
    if ( !_invert ) {
      l_val = _genome.getChromosome(_selected_chromosome).getBlock(_selected_block).pair().len();
      f_val = l_val.floatValue()*sp_start;
      coord = _genome.getChromosome(_selected_chromosome).getBlock(_selected_block).pair().start()+f_val.longValue();
      //coord -= int(_zoom_center-_zoom_sizes[_selected_zoom_level]/2.0);
      //coord /= _zoom_sizes[_selected_zoom_level];
      textAlign( LEFT, BOTTOM );
      text( coord.toString(), ps-text_offset, h );
    }
    else {
      l_val = _genome.getChromosome(_selected_chromosome).getBlock(_selected_block).pair().len();
      f_val = l_val.floatValue()*sp_start;
      coord = _genome.getChromosome(_selected_chromosome).getBlock(_selected_block).pair().stop()+f_val.longValue();
      //coord -= int(_zoom_center-_zoom_sizes[_selected_zoom_level]/2.0);
      //coord /= _zoom_sizes[_selected_zoom_level];
      textAlign( LEFT, BOTTOM );
      text( coord.toString(), ps-text_offset, h );
    }
    
    
    h = _bar_height + 5;
    l_val = _genome.getChromosome(_selected_chromosome).getBlock(_selected_block).len();
    f_val = l_val.floatValue()*sp_start;
    coord = _genome.getChromosome(_selected_chromosome).getBlock(_selected_block).stop()+f_val.longValue();
    //coord -= int(_zoom_center-_zoom_sizes[_selected_zoom_level]/2.0);
    //coord /= _zoom_sizes[_selected_zoom_level];
    textAlign( RIGHT, TOP );
    text( coord.toString(), s+text_offset, h );
    if ( !_invert ) {
      l_val = _genome.getChromosome(_selected_chromosome).getBlock(_selected_block).pair().len();
      f_val = l_val.floatValue()*sp_start;
      coord = _genome.getChromosome(_selected_chromosome).getBlock(_selected_block).pair().stop()+f_val.longValue();
      //coord -= int(_zoom_center-_zoom_sizes[_selected_zoom_level]/2.0);
      //coord /= _zoom_sizes[_selected_zoom_level];
      textAlign( LEFT, TOP );
      text( coord.toString(), ps-text_offset, h );
    }
    else {
      l_val = _genome.getChromosome(_selected_chromosome).getBlock(_selected_block).pair().len();
      f_val = l_val.floatValue()*sp_start;
      coord = _genome.getChromosome(_selected_chromosome).getBlock(_selected_block).pair().start()+f_val.longValue();
      //coord -= int(_zoom_center-_zoom_sizes[_selected_zoom_level]/2.0);
      //coord /= _zoom_sizes[_selected_zoom_level];
      textAlign( LEFT, TOP );
      text( coord.toString(), ps-text_offset, h );
    }
  }
  
  void renderOrientationLegend() {
    float x = -HISTO_BAR_OFFSET-HISTO_BAR_WIDTH;
    float y = _bar_height + _histo_offset + _histo_l + 40;

    //smooth();
    textFont( sans_15 );
    textAlign( LEFT, CENTER );
    strokeWeight( 1.0 );

    //fill( _chrom_color );
    fill( 0 );
    text( "orientation:", x, y );

    x += 20;
    y+= TEXT_DIAMETER/2;

    rectMode( CORNER );
    stroke( _chrom_color );

    fill( _positive_orientation_color );
    rect( x, y, 30, 10 );
    //fill( _chrom_color );
    fill( 0 );
    text( "match", x+40, y );

    y+= TEXT_DIAMETER/2;

    fill( _negative_orientation_color );
    rect( x, y, 30, 10 );
    //fill( _chrom_color );
    fill( 0 );
    text( "inversion", x+40, y );

  }
}


