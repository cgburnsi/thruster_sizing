

if __name__ == '__main__':
    
    props = PropertyTables()
    mech = Mechanism('kinetics.json')
    cfg = ReactorConfig.from_defaults(props)
    cfg.mechanism = mech
    
    result = run_simulation(cfg=cfg, props=props, mech=mech)
    
    plot_summary(result, cfg)
    plot_results(results)
    
    