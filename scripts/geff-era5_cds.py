import cdsapi

c = cdsapi.Client()

c.retrieve(
    'test-fire-pub',
    {
        'format':'zip',
        'model':'reanalysis',
        'variable':'fire_weather_index',
        'year':'2017',
        'month':[
            '01','02','03',
            '04','05','06',
            '07','08','09',
            '10','11','12'
        ],
        'day':[
            '01','02','03',
            '04','05','06',
            '07','08','09',
            '10','11','12',
            '13','14','15',
            '16','17','18',
            '19','20','21',
            '22','23','24',
            '25','26','27',
            '28','29','30',
            '31'
        ]
    },
    'download.zip')