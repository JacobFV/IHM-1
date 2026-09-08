import test from 'node:test';
import assert from 'node:assert/strict';
import {mountTiles} from '../src/tiles.js';

// Minimal DOM surface this widget touches; no browser and no network.
class Element {
 constructor(tag){this.tagName=tag;this.children=[];this.dataset={};this.attributes={};
  this.className='';this.title='';this.disabled=false;this.hidden=false;this._text='';this.innerHTML='';
  this.classList={add:(c)=>{this.className=`${this.className} ${c}`.trim();},
   toggle:(c,on)=>{const has=this.className.split(' ').includes(c);
    if(on??!has)this.classList.add(c);else this.className=this.className.split(' ').filter(x=>x!==c).join(' ');},
   contains:(c)=>this.className.split(' ').includes(c)};}
 set textContent(text){this._text=String(text);this.children=[];}
 get textContent(){return this._text+this.children.map(n=>n.textContent).join('');}
 append(...nodes){this.children.push(...nodes);}
 replaceChildren(...nodes){this._text='';this.children=[...nodes];}
 setAttribute(key,value){this.attributes[key]=String(value);}
 getAttribute(key){return this.attributes[key];}
 click(){this.onclick?.();}
 querySelectorAll(selector){
  const match=(n)=>selector.startsWith('[data-tile=')
   ?n.dataset.tile===selector.slice('[data-tile="'.length,-2)
   :n.tagName===selector;
  return this.children.flatMap(n=>[...(match(n)?[n]:[]),...n.querySelectorAll(selector)]);}
 querySelector(selector){return this.querySelectorAll(selector)[0]||null;}
}

// The catalogue shape the scene tiles actually carry: three exclusive
// environments, and scenes that each name the one they sit on.
const SLOTS=[
 {id:'environment',label:'Environment',exclusive:true,required:true,default:'bed',requires:[]},
 {id:'scene',label:'Scene',exclusive:true,required:false,default:null,requires:[],
  note:'Each scene names the base environment it requires.'},
];
const ITEMS=[
 {id:'studio',label:'Studio',slot:'environment',kind:'environment',requires:[]},
 {id:'floor',label:'Floor',slot:'environment',kind:'environment',requires:[]},
 {id:'bed',label:'Bed',slot:'environment',kind:'environment',requires:[]},
 {id:'bedroom',label:'Bedroom',slot:'scene',kind:'scene',requires:[{slot:'environment',any_of:['bed']}]},
 {id:'grass-field',label:'Grass field',slot:'scene',kind:'scene',requires:[{slot:'environment',any_of:['floor']}]},
 {id:'either-scene',label:'Either',slot:'scene',kind:'scene',requires:[{slot:'environment',any_of:['floor','studio']}]},
];

function host(){
 const changes=[];
 globalThis.document={createElement:(t)=>new Element(t)};
 const node=new Element('div');
 const tiles=mountTiles(node,{onChange:(active,selection)=>changes.push({active,selection})});
 tiles.setItems(ITEMS,{slots:SLOTS,title:'Environment',initial:['bed']});
 return {node,tiles,changes,at:(id)=>node.querySelector(`[data-tile="${id}"]`)};
}

test('a scene on another environment is offered, not disabled, and says what it will also select',()=>{
 const {at}=host();
 assert.equal(at('bedroom').disabled,false,'the scene on the selected environment');
 assert.equal(at('grass-field').disabled,false,'one declared step away, so reachable');
 assert.equal(at('grass-field').title,'Also selects Floor');
 assert.equal(at('bedroom').title,'','nothing to announce when the requirement is already met');
});

test('pressing that scene takes the declared step first, and the environment follows', ()=>{
 const {tiles,changes,at}=host();
 at('grass-field').click();
 assert.deepEqual(tiles.selection,{environment:'floor',scene:'grass-field'});
 assert.deepEqual(changes.at(-1).active.sort(),['floor','grass-field']);
 assert.equal(at('floor').getAttribute('aria-pressed'),'true');
 // And the scene that needed the environment just left comes off with it.
 at('bed').click();
 assert.deepEqual(tiles.selection,{environment:'bed'});
});

test('a requirement with two ways to satisfy it is left to the reader',()=>{
 const {at}=host();
 assert.equal(at('either-scene').disabled,true);
 assert.equal(at('either-scene').title,SLOTS[1].note);
});

test('the reachable route never disturbs the live selection until the press',()=>{
 const {tiles,changes}=host();
 const before={...tiles.selection};
 assert.deepEqual(before,{environment:'bed'});
 assert.deepEqual(tiles.selection,before,'rendering the tiles routed nothing');
 assert.equal(changes.length,0,'and emitted nothing');
});
