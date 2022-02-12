use plotters::prelude::*;
use crate::R;

pub fn plot_1d(x: &[Vec<R>], y: &[Vec<R>], fname: &str,
        title: &str, xlabel: &str, ylabel: &str, labels: &[String],
        col1: (u8, u8, u8), col2: (u8, u8, u8)) 
    -> Result<(), Box<dyn std::error::Error>> 
{
     
    // check that we have the right number of labels
    if labels.len() != y.len() {
        return Err(TextError::new_box(format!("Different numbers of labels and lines ({} and {})",
                            labels.len(), y.len())));
    }
 
    let fname = format!("{}.svg", fname);
    let root = SVGBackend::new(&fname, (640, 480)).into_drawing_area();
        
    // find the minimum and maximum of x and y
    let (x_min, x_max) = find_min_max(x)?;
    let (y_min, y_max) = find_min_max(y)?;
     
    root.fill(&WHITE)?;
    let root = root.margin(20, 20, 20, 20);
    let mut chart = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 20).into_font())
        .x_label_area_size(40)
        .y_label_area_size(50)
        .build_cartesian_2d(x_min..x_max, y_min..y_max)?;

    chart
        .configure_mesh()
        .x_labels(10)
        .y_labels(10)
        .x_desc(xlabel)
        .y_desc(ylabel)
        .axis_desc_style(("sans-serif", 15))
        .draw()?;
    
    let n_plots = y.len();
    for j in 0..n_plots {

        let x = x[j % x.len()].clone();
        let y = y[j].clone();

        // check that x and y have the same size
        if x.len() != y.len() {
            return Err(TextError::new_box(format!("The x and y vectors have different sizes ({} and {})",
                                x.len(), y.len())));
        }
        
        let ep = if n_plots > 1 {
            (j as f64) / (n_plots - 1) as f64
        } else {
            0.
        };
        let color = RGBColor(
            ((col2.0 as f64) * ep + (col1.0 as f64) * (1.-ep)) as u8, 
            ((col2.1 as f64) * ep + (col1.1 as f64) * (1.-ep)) as u8, 
            ((col2.2 as f64) * ep + (col1.2 as f64) * (1.-ep)) as u8,
        );
         
        // build the vector of coordinates
        let mut z = Vec::<(R,R)>::new();
        for i in 0..x.len() {
            z.push((x[i], y[i]));
        }
        let z = z;

        chart.draw_series(LineSeries::new(z.clone(), &color))?
            .label(&labels[j])
            .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &color));
    }

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.5))
        .border_style(&BLACK.mix(0.5))
        .position(plotters::chart::SeriesLabelPosition::LowerRight)
        .draw()?;

    Ok(())
}


fn find_min_max<T: PartialOrd + Copy>(x: &[Vec<T>]) -> Result<(T,T), TextError> {
    
    let mut x_min: Option<T> = None;
    let mut x_max: Option<T> = None;
    
    for j in 0..x.len() {
        for i in 0..x[j].len() {
            if let (Some(x_min_c), Some(x_max_c)) = (x_min, x_max) {
               if x[j][i] < x_min_c { x_min = Some(x[j][i]) };
                if x[j][i] > x_max_c { x_max = Some(x[j][i]) };
            } else {
                x_min = Some(x[j][i]);
                x_max = Some(x[j][i]);
            } 
        }
    }
    
    if let (Some(x_min_c), Some(x_max_c)) = (x_min, x_max) {
        return Ok((x_min_c, x_max_c));
    } else {
        return Err(TextError::new("I can't compute the min and mad of an empty array"
                                  .to_string()));
    }
}


#[derive(Debug, Clone)]
struct TextError { message: String }

impl TextError {
    fn new(s: String) -> TextError {
        TextError { message: s }
    }
    fn new_box(s: String) -> Box<TextError> {
        Box::new(Self::new(s))
    }
}

impl std::fmt::Display for TextError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", &self.message)
    }
}

impl std::error::Error for TextError {} 
