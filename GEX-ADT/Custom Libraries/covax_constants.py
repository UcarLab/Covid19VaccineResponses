import pandas as pd

vax_colors = {
    'J&J' : '#299c0c', #Green 
    'Moderna' : '#1494d9', #Blue
    'Pfizer' : '#d91432' #Red
}

light_vax_colors = {
    'J&J' : '#9bd6f6', #Light blue
    'Moderna' : '#38d510', #Light green
    'Pfizer' : '#f69ba9' #Light red
}

vax_markers = {
    'J&J' : 's',
    'Moderna' : 'o',
    'Pfizer' : 'D'
}



vax_visit_colors = {
    'J&J' : {
        'Baseline 1' : '#006400', # Dark Green
        'Vac1 D1': '#32CD32', # Lime Green
        'Vac1 D7': '#90EE90', # Light Green
    },
    'Moderna' : {
        'Baseline 1': '#000080', # Navy  
        'Vac1 D1': '#0000CD', # Medium Blue  
        'Baseline 2': '#4169E1', # Royal Blue  
        'Vac2 D1': '#1E90FF', # Dodger Blue  
        'Vac2 D7': '#87CEEB', # Sky Blue  
        
    },
    'Pfizer' : {

        'Baseline 1': '#8B0000', # Dark Red  
        'Vac1 D1': '#B22222', # Firebrick  
        'Baseline 2': '#DC143C', # Crimson  
        'Vac2 D1': '#CD5C5C', # Indian Red
        'Vac2 D7': '#FF6347', # Tomato  
        
    }
}
