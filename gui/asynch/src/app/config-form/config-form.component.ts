import { Input, Component, OnInit } from '@angular/core';
import { ForcingStateConfig, ForcingConfig, AsynchConfig } from '../config';
import { ModelMeta,  models } from '../models';


@Component({
  selector: 'app-config-form',
  templateUrl: './config-form.component.html',
  styleUrls: ['./config-form.component.scss']
})
export class ConfigFormComponent implements OnInit {

  models: ModelMeta[] = models;
  outputFunctions: string[] = [
    'Time', 'TimeI', 'State0', 'State1', 'State2', 'State3', 'State4'
  ]
  peakFunctions: string[] = [
    'Classic', 'Forecast'
  ]

  private _config: AsynchConfig;  
  private _selectedModel: ModelMeta;
  
  constructor() {
    //this._config = new AsynchConfig();
  }

  ngOnInit() {
  }
  
  @Input()
  set selectedModel(model: ModelMeta) {
    if (model) {
      this._selectedModel = model;
      this._config.init(model);
    }
  }
  
  get selectedModel(): ModelMeta { return this._selectedModel; }
  
  @Input()
  set config(config: AsynchConfig) {
    if (config) {
      this._config = config;
      this._selectedModel = models.find((model) => {return model.uid === config.model;})
    } else {
      this._selectedModel = null;
    }
  }
  
  get config(): AsynchConfig { return this._config; }
  
  setFilename(file: File, filename: string) {
    filename = file.name;
  }
  
  toggleForcing(i: number) {
    if (this._config.forcings.timeseries[i]) {
      this._config.forcings.timeseries[i] = null;
    } else {
      this._config.forcings.timeseries[i] = new ForcingConfig();
    }
  }
  
  toggleState() {
    if (this._config.forcings.state) {
      this._config.forcings.state = null;
    } else {
      this._config.forcings.state = new ForcingStateConfig();
    }
  }
  
  useDefaultSolverConfig(useDefault: boolean) {
    if (useDefault) {
      this.config.removeSolver();
    } else  {
      this.config.addSolver();
    }
      
  }
}
