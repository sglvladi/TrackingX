lambda = 5;
limits = [0 10; 0 10];

clutter_model = PoissonRateUniformPositionX('ClutterRate',lambda,'Limits',limits);

lambda_samples = clutter_model.random(1000,'cardinality');
spatial_samples = clutter_model.random();