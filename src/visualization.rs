use plotters::prelude::*;

pub fn plot(
    path: Vec<(Vec<usize>, Vec<usize>)>,
    len_a: usize,
    len_b: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    // Make this into a dataframe
    // let mut x: Vec<i32> = Vec::new();
    // let mut y: Vec<i32> = Vec::new();
    let mut x_y: Vec<(f32, f32)> = Vec::new();
    for (afp_a, afp_b) in &path {
        for (_x, _y) in afp_a.iter().zip(afp_b.iter()) {
            // x.push(*_x as i32);
            // y.push(*_y as i32);
            // for i in 1..8 {
            //     x_y.push((*_x as f32 + i as f32, *_y as f32 + i as f32));
            // }
            x_y.push((*_x as f32, *_y as f32));
        }

        // x_y.push((i as f32, j as f32));

        // for _x in i..afp_a[0] + window_size {
        //     x.push(_x as i32);
        // }
        // for _y in j..afp_b[0] + window_size {
        //     y.push(_y as i32);
        // }
    }

    // let df = DataFrame::new(vec![Series::new("x", x), Series::new("y", y)]).unwrap();
    // println!("{:?}", df);

    // println!("{:?}", best_path);

    let root = BitMapBackend::new("debug/plot.png", (640, 480)).into_drawing_area();
    let _ = root.fill(&WHITE);
    let root = root.margin(10, 10, 10, 10);
    // After this point, we should be able to construct a chart context
    let mut chart = ChartBuilder::on(&root)
        // Set the caption of the chart
        // .caption("This is our first plot", ("sans-serif", 40).into_font())
        // Set the size of the label region
        .x_label_area_size(20)
        .y_label_area_size(20)
        // Finally attach a coordinate on the drawing area and make a chart context
        .build_cartesian_2d(0f32..len_a as f32, 0f32..len_b as f32)?;

    // Then we can draw a mesh
    chart
        .configure_mesh()
        // We can customize the maximum number of labels allowed for each axis
        .x_labels(10)
        .y_labels(10)
        // We can also change the format of the label text
        .y_label_formatter(&|x| format!("{:.0}", x))
        .x_label_formatter(&|x| format!("{:.0}", x))
        .draw()?;

    // Sort x_y by the first element
    x_y.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    // for (x, y) in &x_y {
    //     println!("x: {}, y: {}", x, y);
    // }

    // And we can draw something in the drawing area
    chart.draw_series(LineSeries::new(x_y.clone(), &BLACK))?;

    // Similarly, we can draw point series
    chart.draw_series(PointSeries::of_element(x_y, 1, &BLACK, &|c, s, st| {
        EmptyElement::at(c)    // We want to construct a composed element on-the-fly
          + Circle::new((0,0),s,st.filled() )
    }))?;
    root.present()?;

    Ok(())
}
