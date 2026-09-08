import {chromium} from 'playwright';
import {mkdirSync,writeFileSync} from 'node:fs';
import {spawn} from 'node:child_process';
const server=spawn('.venv/bin/python',['-m','ihm','serve','--port','8766'],{env:{...process.env,PYTHONPATH:'.'},stdio:['ignore','pipe','inherit']});
await new Promise((resolve,reject)=>{server.stdout.once('data',resolve);server.once('exit',code=>reject(Error('Server exited '+code)));});
process.on('exit',()=>server.kill());
const out='test-results/environment-upgrade';mkdirSync(out,{recursive:true});
const browser=await chromium.launch({executablePath:'/usr/bin/google-chrome',headless:true,args:['--no-sandbox','--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader']});
const page=await browser.newPage({viewport:{width:1600,height:1050}});const errors=[];
page.on('pageerror',e=>errors.push(e.message));page.on('response',r=>{if(r.status()>=400&&/api\/scene/.test(r.url()))errors.push(r.status()+' '+r.url());});
await page.goto('http://127.0.0.1:8766');
await page.locator('#environment-tiles [data-tile="bedroom"]').waitFor({timeout:120000});
for(const id of (process.env.ENVIRONMENT_SCENES?.split(',') || ['bedroom','hospital-room','clinic-room','play-floor','grass-field','garden-patio'])) {
 if(['clinic-room','play-floor','grass-field','garden-patio'].includes(id))await page.locator('#environment-tiles [data-tile="floor"]').click();
 await page.locator(`#environment-tiles [data-tile="${id}"]`).click();
 await page.waitForTimeout(1800);
 await page.locator('#scene').screenshot({path:`${out}/${id}.png`});
}
writeFileSync(`${out}/browser-errors.json`,JSON.stringify(errors,null,2));
await browser.close();server.kill();if(errors.length)throw Error(errors.join('\n'));
